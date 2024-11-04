// #include "main.h"

#include <iostream>
#include <string>

#include "base_cpp/output.h"
#include <molecule/molecule.h>
#include <base_cpp/scanner.h>
#include <molecule/canonical_smiles_saver.h>
#include <molecule/molecule_auto_loader.h>


using namespace indigo;

void loadMolecule(const char* buf, Molecule& molecule)
{
    BufferScanner scanner(buf);
    MoleculeAutoLoader loader(scanner);
    loader.loadMolecule(molecule);
}

std::string smiles(Molecule& m)
{
    Array<char> smiles;
    ArrayOutput output(smiles);
    SmilesSaver saver(output);
    saver.saveMolecule(m);
    return {smiles.ptr(), static_cast<std::size_t>(smiles.size())};
}


int main() {
    Molecule mol;
    Array<char> descr;
    loadMolecule("C1OCC1", mol);
    for (int v : mol.vertices()) {
        int atom_num = mol.getAtomNumber(v);
        mol.getAtomDescription(v, descr);
        std::cout << atom_num << descr << std::endl;
    }

    // loadMolecule(CAFFEINE, mol);
    std::cout << smiles(mol) << std::endl;
    return 0;
}
