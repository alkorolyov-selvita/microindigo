//
// Created by ergot on 11/3/24.
//

#ifndef MAIN_H
#define MAIN_H

#endif //MAIN_H

#define DATA_PATH "/home/ergot/projects/microindigo/data"

#include <string>

#include <molecule/molecule.h>
#include <reaction/reaction.h>
#include <reaction/query_reaction.h>

#define METHANE "C"
#define BENZENE "C1=CC=CC=C1"
#define BENZENE_AROMATIC "c1ccccc1"
#define CAFFEINE "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
#define SULFASALAZINE "C1=CC=NC(=C1)NS(=O)(=O)C2=CC=C(C=C2)N=NC3=CC(=C(C=C3)O)C(=O)O"

namespace indigo
{
    static void loadMolecule(const char* buf, Molecule& molecule);
    static std::string smiles(Molecule& m);
}

