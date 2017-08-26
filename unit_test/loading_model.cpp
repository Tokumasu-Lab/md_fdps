//***************************************************************************************
//  This is unit test of loading model parameter.
//***************************************************************************************

#include <tuple>
#include <iostream>

#include <particle_simulator.hpp>

#include "md_fdps_enum_model.hpp"
#include "md_fdps_coef_table.hpp"
#include "md_fdps_atom_class.hpp"
#include "md_fdps_loading_condition.hpp"
#include "md_fdps_loading_model.hpp"

int main(int argc, char *argv[]) {

    PS::Initialize(argc, argv);

    if(PS::Comm::getRank() == 0) {
        //--- molecular model loading test
        std::cout << "TEST: molecular model loading test..." << std::endl;

        if(argc <= 1){
            std::cout << "    no input." << std::endl
                      << "    usage is:  $ ./test_model.x [model_file_name]";
        }

        System::model_list.clear();
        for(int i=1; i<argc; i++){

            std::string model_name = argv[i];
            std::cout << std::endl;
            std::cout << "  model file: " << model_name << std::endl;

            MolName model = ENUM::which_MolName(model_name);
            System::model_list.push_back( std::make_pair(model, 0) );
            System::model_template.push_back( std::vector<Atom_FP>{} );
            System::bond_template.push_back( std::vector<std::vector<PS::S64>>{} );

            System::model_template.at(i-1).clear();
            System::bond_template.at(i-1).clear();

            loading_model_parameter(model_name,
                                    System::model_template.at(i-1),
                                    System::bond_template.at(i-1),
                                    MODEL::coefTable_elem,
                                    MODEL::coefTable_bond,
                                    MODEL::coefTable_angle,
                                    MODEL::coefTable_torsion);

            //--- show result
            std::cout << endl;
            std::cout << "result: atom" << std::endl;
            std::cout << "  \"MolID\" was set as illigal value(-1). this is model template." << std::endl;
            MODEL::print_model_template( System::model_template.at(i-1),
                                         System::bond_template.at(i-1) );


            std::cout << "result: inter" << std::endl;
            MODEL::print_coef_table( MODEL::coefTable_elem,
                                     ENUM::which_MolName(model_name) );


            std::cout << "result: bond" << std::endl;
            MODEL::print_coef_table( MODEL::coefTable_bond,
                                     ENUM::which_MolName(model_name) );


            std::cout << "result: angle" << std::endl;
            MODEL::print_coef_table( MODEL::coefTable_angle,
                                     ENUM::which_MolName(model_name) );


            std::cout << "result: torsion" << std::endl;
            MODEL::print_coef_table( MODEL::coefTable_torsion,
                                     ENUM::which_MolName(model_name) );
        }

        std::cout << "    the test succeeded." << std::endl;
    }

    PS::Finalize();
    return 0;
}
