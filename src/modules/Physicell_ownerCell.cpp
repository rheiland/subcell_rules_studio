/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include <algorithm>
#include "PhysiCell_ownerCell.h"
#include "../core/PhysiCell_phenotype.h"
#include <csignal>

namespace PhysiCell{

double persistence_time = 100; // TODO: set this in XML

//owner cells are created with an id, list of member agents, target volume
// and growth rate, and cell definition with phenotype
Owner_cell::Owner_cell(int _id, Cell_Definition* cd)
{
	id = _id; 
    deathVolume = parameters.doubles("death_volume");
    targetVolume = cd->phenotype.volume.total;
    volumetric_growth_rate = 0;
    cell_definition = cd;
    std::cout << "Constructor called " << id << std::endl; 
};

//master list of owner cells in simulation
std::vector<Owner_cell*> owner_cells;


int get_max_owner_id()
{
    //give the agent a unique ID  
    static int max_owner_cell_ID = 0; //gets set to 0 the first time, but keeps track persistently; not thread-safe
    auto id = max_owner_cell_ID;
    max_owner_cell_ID++; 
    return id + max_init_type_id;
}

Cell_Definition* initialize_owner_cell_definition( int id, std::string parent_cell_definition_name )
{
	return initialize_cell_definition(
		parent_cell_definition_name + " " + std::to_string(id), 
		3 + id, 
		parent_cell_definition_name,
		pugi::xml_node(),
		pugi::xml_node()
	);
}

Cell_Definition* initialize_member_cell_definition( int id, std::string parent_cell_definition_name )
{
	Cell_Definition* result = initialize_cell_definition(
		parent_cell_definition_name + " " + std::to_string(3 + id + MAX_OWNER_CELLS), 
		3 + id + MAX_OWNER_CELLS, 
		parent_cell_definition_name,
		pugi::xml_node(),
		pugi::xml_node()
	);
    result->phenotype.volume.total = subcell_volume;
    return result;
}

void initialize_subcell_definitions_from_pugixml( pugi::xml_node root )
{
	// Get CSV path and parent cell type name from XML
	std::tuple<std::string, std::string, std::string, double> init_conditions_data = parse_initial_conditions_from_pugixml(root);
    owner_cell_definition_name = std::get<1>(init_conditions_data);
    member_cell_definition_name = std::get<2>(init_conditions_data);
    double subcell_radius = std::get<3>(init_conditions_data);
    subcell_volume = 4 / 3 * 3.14 * pow(subcell_radius, 3.0);
	std::cout << "Creating cell definitions from CSV file " << std::get<0>(init_conditions_data) << " ... " << std::endl; 

	std::ifstream file( std::get<0>(init_conditions_data), std::ios::in );
	if( !file )
	{ 
		std::cout << "Error: " << std::get<0>(init_conditions_data) << " not found during cell loading. Quitting." << std::endl; 
		exit(-1);
	}

	std::string line;
	std::vector<int> type_ids;
	while (std::getline(file, line))
	{
		std::vector<double> data;
		csv_to_vector( line.c_str() , data ); 

		if( data.size() != 4 )
		{
			std::cout << "Error! Importing cells from a CSV file expects each row to be x,y,z,typeID." << std::endl;
			exit(-1);
		}

		int type_id = (int) data[3]; 
		if (std::find(type_ids.begin(), type_ids.end(), type_id) == type_ids.end())
		{
            initialize_owner_cell_definition(type_id, std::get<1>(init_conditions_data));
			initialize_member_cell_definition(type_id, std::get<2>(init_conditions_data));
			type_ids.push_back(type_id);
		}
	}

    int max_type_id = 0;
    for (auto it : type_ids) 
    {
        if (it > max_type_id)
        {
            max_type_id = it;
        }
    }
    max_init_type_id = max_type_id;

	file.close(); 

}

void initialize_subcell_definitions_from_pugixml( void )
{
	initialize_subcell_definitions_from_pugixml( physicell_config_root );
	return; 
}


Owner_cell* find_ownercell( int search_id, Cell_Definition* cd)
{
	Owner_cell* output = NULL;

    // loop through owner_cells and break the loop if id is in this list
	// #pragma omp parallel for
    // possible solution for above: https://stackoverflow.com/questions/9793791/parallel-openmp-loop-with-break-statement
    for( int n = 0 ; n < owner_cells.size() ; n++ )
    {
        if(owner_cells[n]->id == search_id)
        {
            output = owner_cells[n];
        }
    }
	
    // if owner cell with this id doesn't exsit yet, create one
	if( output == NULL )
	{
        std::cout << "Creating new owner cell: " << search_id << std::endl; 
		output = new Owner_cell( search_id, cd);	
        owner_cells.push_back( output );
	}
	
	return output; 	
}

void handle_member_cells_division()
{
    // add any new daughter members to their mother member's owner cell
    for (auto it : cells_that_just_divided) {
        auto owner = owner_cell_for_member[std::get<1>(it)];
        if (!owner)
        {
            std::cout << "Error! Failed to find owner cell for member[" << std::get<1>(it)->ID << "]" << std::endl; 
            exit(-1);
        }
        owner->add_member( std::get<2>(it) );
    }
}

void Owner_cell::remove_dead_members()
{
    // remove any dead member cells from their owner cell
    for( int n = 0 ; n < member_cells.size() ; n++ )
    { 
        if (member_cells[n]->phenotype.death.dead)
        {
            remove_member( member_cells[n] );
        }
    }
}

//timestep update for owner cells
void update_owner_cells( double dt )
{
    handle_member_cells_division();
    
    for( int n = 0 ; n < owner_cells.size() ; n++ )
    { 
        owner_cells[n]->update( dt ); 
    }
    return;
}

//add member agent to list of agents belonging to owner cell
void Owner_cell::add_member( Cell* member_cell )
{
    member_cells.push_back( member_cell );
    owner_cell_for_member[member_cell] = this;
}

//remove the member agent from the list of agents belonging to owner cell
void Owner_cell::remove_member( Cell* member_cell )
{
    auto it = std::find(member_cells.begin(), member_cells.end(), member_cell );
    if(it != member_cells.end())
        member_cells.erase(it);
    else{
        std::cout << "Error: Owner " << id << " failed to remove member[" << member_cell->ID << "]" << std::endl; 
        exit(-1);
    }
    owner_cell_for_member.erase( member_cell );
}

double Owner_cell::get_volume()
{
    // add up subcell volumes to get current total cell volume and return it
    double total_volume = 0;
    for( int n = 0 ; n < member_cells.size() ; n++ )
    {
        total_volume += member_cells[n]->phenotype.volume.total;
    }
    return total_volume;
}

// get distance between two 3D points
double distance(std::vector<double> position1, std::vector<double> position2)
{
    return sqrt(pow(position2[0] - position1[0], 2) +
                pow(position2[1] - position1[1], 2) +
                pow(position2[2] - position1[2], 2) * 1.0);
}

//calculate dot product of two 3D vectors
double get_dot_product(std::vector<double> a, std::vector<double> b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

//calculate difference between two 3D vectors
std::vector<double> subtract_vectors(std::vector<double> a, std::vector<double> b)
{
    try {
        std::cout << "a = (" << a[0] << ", " << a[1] << ", " << a[2] << ")" << std::endl; 
    }
    catch (const std::exception&) {
        std::cout << "Vector A has a problem" << std::endl; 
    }

    try {
        std::cout << "b = (" << b[0] << ", " << b[1] << ", " << b[2] << ")" << std::endl; 
    }
    catch (const std::exception&) {
        std::cout << "Vector B has a problem" << std::endl; 
    }

    try {
        std::vector<double> x = { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
    }
    catch (const std::exception&) {
        std::cout << "subtraction failed!!!" << std::endl; 
    }

    return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
}

//calculate product of two 3D vectors
std::vector<double> multiply_vector(std::vector<double> v, double factor)
{
    return { v[0] * factor, v[1] * factor, v[2] * factor };
}

//calculate difference between two 3D vectors
std::vector<double> divide_vector(std::vector<double> v, double divisor)
{
    return { v[0] / divisor, v[1] / divisor, v[2] / divisor };
}

//find plane in 3D along which the owner cell will divide
std::vector<double> get_dividing_plane_point(std::vector<double> long_axis_normal, std::vector<double> centroid, std::vector<double> position1)
{
    auto v = subtract_vectors(centroid, position1);
    return subtract_vectors(centroid, subtract_vectors(v, multiply_vector(long_axis_normal, get_dot_product(v, long_axis_normal))));
}

// determine whether a member cell will be part of the new daughter cell
bool test_if_subcell_moves_to_daughter( 
    std::vector<Cell*> member_cells, 
    int member_index, 
    std::vector<double> long_axis_normal, 
    std::vector<double> dividing_plane_point
)
{
    std::vector<double> test_point = member_cells[member_index]->position;
    double plane_side = get_dot_product(long_axis_normal, test_point) - get_dot_product(long_axis_normal, dividing_plane_point);
    return plane_side > 0;
}

// split up the subcells into two daughter ownercells
std::vector<Cell*> get_subcells_to_move_to_daughter( std::vector<Cell*> member_cells )
{
    // get list of member positions
    std::vector<std::vector<double>> positions;
    positions.clear(); //make sure this is empty first
    for( int n = 0 ; n < member_cells.size() ; n++ )
    { 
        positions.push_back( member_cells[n]->position );
    }

    // find long axis of the positions cloud
    double max_distance;
    std::vector<double> position1;
    std::vector<double> position2;
    //#pragma omp parallel for
    for( int n = 0 ; n < positions.size() ; n++ )
    { 
        for( int m = n + 1 ; m < positions.size() ; m++ )
        { 
            // TODO: replace with squared distances, fnct exists for this already in BioFVM vector
            // TODO: also look there for replacements for vector math helper fncts
            double d = distance(positions[n], positions[m]);
            if (d > max_distance) 
            {
                //#pragma omp critical
                {
                    max_distance = d;
                    position1 = positions[n];
                    position2 = positions[m];
                }
            }
        }
    }
    std::vector<double> long_axis = subtract_vectors(position2, position1);

    // find centroid point of the positions cloud
    std::vector<double> centroid = {0, 0, 0};
    for( int n = 0 ; n < positions.size() ; n++ )
    {
        centroid[0] += positions[n][0];
        centroid[1] += positions[n][1];
        centroid[2] += positions[n][2];
    }
    centroid /= positions.size();
    
    // get plane perpendicular to axis through centroid
    // ref: https://web.ma.utexas.edu/users/m408m/Display12-5-3.shtml
    std::vector<double> long_axis_normal = divide_vector(long_axis, sqrt(get_dot_product(long_axis, long_axis)));
    std::vector<double> plane_point = get_dividing_plane_point(long_axis_normal, centroid, position1);

    // sort subcells on either side of plane
    std::vector<Cell*> cells_to_move;
    cells_to_move.clear();
    //#pragma omp parallel for
    for( int n = 0 ; n < member_cells.size() ; n++ )
    {
        if (test_if_subcell_moves_to_daughter(member_cells, n, long_axis_normal, plane_point))
        {
            //#pragma omp critical
            {
                cells_to_move.push_back(member_cells[n]);
            }
        }
    }
    return cells_to_move;
}

void Owner_cell::die( )
{
    // tell member cells to die: start death and remove from owner cell map
    for( int n = 0 ; n < member_cells.size() ; n++ )
    {
        member_cells[n]->start_death(0);
        owner_cell_for_member.erase( member_cells[n] );
    }

    // external list of owner cells needs to remove this owner cell from its list
    auto it = std::find(owner_cells.begin(), owner_cells.end(), this );
    if(it != owner_cells.end())
        owner_cells.erase(it);
    else{
        std::cout << "Error: Owner " << id << " failed to die" << std::endl; 
        exit(-1);
    }

    //delete this owner cell object
    delete this;
}


Owner_cell* Owner_cell::divide( )
{ 
    // create new owner cell and move half of subcells to it
    std::vector<Cell*> cells_to_move = get_subcells_to_move_to_daughter( member_cells );
    auto newOwnerType = get_max_owner_id();
    Cell_Definition* newOwnerCD = initialize_owner_cell_definition( 
        newOwnerType, owner_cell_definition_name );
    Cell_Definition* newMemberCD = initialize_member_cell_definition( 
        newOwnerType + MAX_OWNER_CELLS, member_cell_definition_name );
    auto newOwner = new Owner_cell( newOwnerType, newOwnerCD );
    for( int n = 0 ; n < cells_to_move.size() ; n++ )
    {
        remove_member( cells_to_move[n] );
        cells_to_move[n]->convert_to_cell_definition( *newMemberCD );
        newOwner->add_member( cells_to_move[n] );
    }
    newOwner->initialize_motility();
    return newOwner;
} 

std::vector<double> Owner_cell::get_random_vector(bool restrict_to_2D)
{
    // choose a uniformly random unit vector 
    double temp_angle = 6.28318530717959*UniformRandom();
    double temp_phi = 3.1415926535897932384626433832795*UniformRandom();
    
    double sin_phi = sin(temp_phi);
    double cos_phi = cos(temp_phi);
    
    if( restrict_to_2D == true )
    {
        sin_phi = 1.0; 
        cos_phi = 0.0;
    }
    
    std::vector<double> randvec; 
    randvec.resize(3,sin_phi); 
    
    randvec[0] *= cos( temp_angle ); // cos(theta)*sin(phi)
    randvec[1] *= sin( temp_angle ); // sin(theta)*sin(phi) 
    randvec[2] = cos_phi; //  cos(phi)
    return randvec;
}

void Owner_cell::initialize_motility()
{
    //check if motile and if not don't do anything
    if( cell_definition->phenotype.motility.is_motile == false )
	{
		cell_definition->phenotype.motility.motility_vector.assign( 3, 0.0 ); 
		return; 
	}

    std::cout << "Motility True - Initializing motility " << id << std::endl; 

    //choose random direction to initialize subcell migration direction
    std::vector<double> randvec = get_random_vector(true);

    //initialize motility coordianted across subcells for persistence (0-1) or persistence time?
    for( int n = 0 ; n < member_cells.size() ; n++ )
    {
        member_cells[n]->phenotype.motility.migration_bias_direction = randvec;
        cell_definition->phenotype.motility.migration_bias = 1;
    }
    last_direction_change_time = PhysiCell_globals.current_time;
}


void Owner_cell::update_motility()
{
     //check if motile and if not don't do anything
    if( cell_definition->phenotype.motility.is_motile == false )
	{
		cell_definition->phenotype.motility.motility_vector.assign( 3, 0.0 ); 
		return; 
	}
    // check if enough time has passed to change direction
    if( PhysiCell_globals.current_time - last_direction_change_time < persistence_time )
    {
        return;
    }
    last_direction_change_time = PhysiCell_globals.current_time;
    //updates coordinate subcells to have same motility vector
    std::vector<double> randvec = get_random_vector(true);
    for( int n = 0 ; n < member_cells.size() ; n++ )
    {
        member_cells[n]->phenotype.motility.migration_bias_direction = randvec;
        cell_definition->phenotype.motility.migration_bias = 1;
    }
}

void set_initial_member_parameters()
{
    std::cout << "Initialize member parameters for all!!!!!!!!!!!!!!!!" << std::endl; 
    for( int n = 0 ; n < owner_cells.size() ; n++ )
    {
        owner_cells[n]->initialize_motility();
    }
}

// just update division, since we're using the existing volume control machinery
void Owner_cell::update( double dt )
{
    remove_dead_members();

    // get current volume to check against death and division
    auto current_volume = get_volume();

    // is it time to die?
    if ( current_volume < deathVolume )
    {
        std::cout << "DIE owner[" << id << "] " << std::endl; 
        die();

    }

    // is it time to divide?
    else if( current_volume > targetVolume )
    {
        auto daughter_cell = divide();
        owner_cells.push_back( daughter_cell );
        std::cout << "DIVIDE owner[" << id << "] -> owner[" << id << "] + owner[" << daughter_cell->id << "]" << std::endl; 
    }

    update_motility();
} 

};