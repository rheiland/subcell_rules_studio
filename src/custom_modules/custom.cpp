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

#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	
	cell_defaults.functions.volume_update_function = NULL; // standard_volume_update_function;
	cell_defaults.functions.update_velocity = heterotypic_update_cell_velocity; // standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   Create a cell definition for each cell ID in initial conditions CSV 
	   specified in the XML config file. 
	*/
	initialize_subcell_definitions_from_pugixml();
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
		
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void arrange_cells_random( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// add some margin 
	
	double margin = parameters.doubles("initial_margin"); 
	
	Xmin += margin * Xrange; 
	Xmax -= margin * Xrange;
	Xrange *= (1-2*margin); 

	Ymin += margin * Yrange; 
	Ymax -= margin * Yrange;
	Yrange *= (1-2*margin); 
	
	Zmin += margin * Zrange; 
	Zmax -= margin * Zrange;
	Zrange *= (1-2*margin); 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	// place default cells 
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
}

void arrange_cell_blobs( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// add some margin 
	
	double margin = parameters.doubles("initial_margin"); 
	
	Xmin += margin * Xrange; 
	Xmax -= margin * Xrange;
	Xrange *= (1-2*margin); 

	Ymin += margin * Yrange; 
	Ymax -= margin * Yrange;
	Yrange *= (1-2*margin); 
	
	Zmin += margin * Zrange; 
	Zmax -= margin * Zrange;
	Zrange *= (1-2*margin); 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	// place default cells 
	
	int number_of_types = cell_definitions_by_index.size();
	
	double x = Xmin; 
	double y = Ymin; 
	double spacing = 120; 
	double offset = 0.5*spacing; 
	
	int cell_type = 0; 
	
	int number_of_cells_per_block = 60; // 100 
	int row_number = 0; 
	
	while( y < Ymax )
	{
		while( x < Xmax )
		{
			// select type 
			Cell_Definition* pCD = cell_definitions_by_index[cell_type]; 
			
			// place cells
			
			for( int n=0 ; n < number_of_cells_per_block ; n++ )
			{
				std::vector<double> position = { 0,0,0 }; 
				position[0] = x + UniformRandom()*spacing; 
				position[1] = y + UniformRandom()*spacing; 
				Cell* pC = create_cell( *pCD ); 
				pC->assign_position( position );
			}
			
			// move start position 
			x += (0.8*spacing); 
			
			// increment cell type 
			cell_type++; 
			if( cell_type >= number_of_types )
			{ cell_type = 0; } 
		}
		row_number++; 
		x = Xmin; 
		if( row_number % 2 == 1 )
		{ x += offset; }
		
		y += (0.8*spacing); 
	}
}


void setup_tissue( void )
{
	// arrange_cell_blobs(); 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring_wrapped(pCell); } 
// return paint_by_number_cell_coloring(pCell); }



void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 


void add_heterotypic_potentials(Cell* my_cell , Cell* other_agent)
{
	// if( this->ID == other_agent->ID )
	if( my_cell == other_agent )
	{ return; }
	
	static double rel_heterotypic_adhesion = parameters.doubles("relative_heterotypic_adhesion"); 

	// 12 uniform neighbors at a close packing distance, after dividing out all constants
	static double simple_pressure_scale = 0.027288820670331; // 12 * (1 - sqrt(pi/(2*sqrt(3))))^2 
	// 9.820170012151277; // 12 * ( 1 - sqrt(2*pi/sqrt(3)))^2

	double distance = 0; 
	for( int i = 0 ; i < 3 ; i++ ) 
	{ 
		my_cell->displacement[i] = my_cell->position[i] - (*other_agent).position[i]; 
		distance += (my_cell->displacement[i]) * (my_cell->displacement[i]); 
	}
	// Make sure that the distance is not zero
	
	distance = std::max(sqrt(distance), 0.00001); 
	
	//Repulsive
	double R = my_cell->phenotype.geometry.radius + (*other_agent).phenotype.geometry.radius; 
	
	double RN = my_cell->phenotype.geometry.nuclear_radius + (*other_agent).phenotype.geometry.nuclear_radius;	
	double temp_r, c;
	if( distance > R ) 
	{
		temp_r=0;
	}
	else
	{
		temp_r = -distance; // -d
		temp_r /= R; // -d/R
		temp_r += 1.0; // 1-d/R
		temp_r *= temp_r; // (1-d/R)^2 
		
		// add the relative pressure contribution 
		my_cell->state.simple_pressure += ( temp_r / simple_pressure_scale ); // New July 2017 
	}
	
	// August 2017 - back to the original if both have same coefficient 
	double effective_repulsion = sqrt( my_cell->phenotype.mechanics.cell_cell_repulsion_strength * other_agent->phenotype.mechanics.cell_cell_repulsion_strength );
	temp_r *= effective_repulsion; 
	
	//////////////////////////////////////////////////////////////////
	
	// Adhesive
	double max_interactive_distance = my_cell->phenotype.mechanics.relative_maximum_adhesion_distance * my_cell->phenotype.geometry.radius + 
		(*other_agent).phenotype.mechanics.relative_maximum_adhesion_distance * (*other_agent).phenotype.geometry.radius;
		
	if(distance < max_interactive_distance ) 
	{	
		// double temp_a = 1 - distance/max_interactive_distance; 
		double temp_a = -distance; // -d
		temp_a /= max_interactive_distance; // -d/S
		temp_a += 1.0; // 1 - d/S 
		temp_a *= temp_a; // (1-d/S)^2 
		
		// August 2017 - back to the original if both have same coefficient 
		double effective_adhesion = sqrt( my_cell->phenotype.mechanics.cell_cell_adhesion_strength * other_agent->phenotype.mechanics.cell_cell_adhesion_strength ); 
		if( my_cell->type != other_agent->type )
		{ effective_adhesion *= rel_heterotypic_adhesion; }
		temp_a *= effective_adhesion; 
		
		temp_r -= temp_a;
	}
	/////////////////////////////////////////////////////////////////
	if( fabs(temp_r) < 1e-16 )
	{ return; }
	temp_r /= distance;
	axpy( &(my_cell->velocity) , temp_r , my_cell->displacement ); 
	
	return;
}

void heterotypic_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
	if( pCell->functions.add_cell_basement_membrane_interactions )
	{
		pCell->functions.add_cell_basement_membrane_interactions(pCell, phenotype,dt);
	}
	
	pCell->state.simple_pressure = 0.0; 
	
	//First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		// pCell->add_potentials(*neighbor);
		add_heterotypic_potentials( pCell, *neighbor ); 
	}
	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{
			// pCell->add_potentials(*neighbor);
			add_heterotypic_potentials( pCell, *neighbor ); 
		}
	}

	pCell->update_motility_vector(dt); 
	pCell->velocity += phenotype.motility.motility_vector; 
	
	return; 
}


std::vector<std::string> paint_by_number_cell_coloring_wrapped( Cell* pCell )
{
	static std::vector< std::string > colors(0); 
	static bool setup_done = false; 
	if( setup_done == false )
	{
		colors.push_back( "grey" ); // default color will be grey 

		colors.push_back( "red" );
		colors.push_back( "yellow" ); 	
		colors.push_back( "green" ); 	
		colors.push_back( "blue" ); // 4
		
		colors.push_back( "magenta" ); 	
		colors.push_back( "orange" ); 	
		colors.push_back( "lime" ); 	
		colors.push_back( "cyan" ); // 8
		
		colors.push_back( "hotpink" ); 	
		colors.push_back( "peachpuff" ); 	
		colors.push_back( "darkseagreen" ); 	
		colors.push_back( "lightskyblue" ); // 12

		colors.push_back( "darkred" );
		colors.push_back( "goldenrod" ); 	
		colors.push_back( "darkgreen" ); 	
		colors.push_back( "darkblue" );  // 16
		
		setup_done = true; 
	}
	
	int number_of_colors = colors.size(); 
	
	// start all black 
	
	std::vector<std::string> output = { "black", "black", "black", "black" }; 
	
	// paint by number -- by cell type 
	
	// modular arithmetic  
	// int type = (int) ( ((unsigned int) pCell->type) % 13 ); 
	int agent_type = pCell->type % number_of_colors; 
	
//	std::cout << type << std::endl; 
	
	std::string interior_color = "white"; 
	if( agent_type < number_of_colors )
	{ interior_color = colors[ agent_type ]; }
	
	output[0] = interior_color; // set cytoplasm color 
	
	if( pCell->phenotype.death.dead == false ) // if live, color nucleus same color 
	{
		output[2] = interior_color; 
		output[3] = interior_color; 
	}
	else
	{
		// apoptotic cells will retain a black nucleus 
		// if necrotic, color the nucleus brown 
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
		{
			output[2] = "rgb(139,69,19)";
			output[3] = "rgb(139,69,19)";
		}
	}
	
	return output; 
}	
