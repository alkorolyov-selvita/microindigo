// #include "main.h"

#include "main.h"

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

constexpr int ONEHOT_SIZE = 10;
void onehot_encode(int atom_num, std::array<int, ONEHOT_SIZE>& onehot) {
    int onehot_idx[256] = {9}; // fill with 9 for unknown atom first
    // define known atoms
    onehot_idx[6] = 0;  // C
    onehot_idx[7] = 1;  // N
    onehot_idx[8] = 2;  // O
    onehot_idx[15] = 3; // P
    onehot_idx[16] = 4; // S
    onehot_idx[9] = 5;  // F
    onehot_idx[17] = 6; // Cl
    onehot_idx[35] = 7; // Br
    onehot_idx[53] = 8; // I
    int idx = onehot_idx[atom_num];

    onehot.fill(0);
    onehot[idx] = 1;
}

template <std::size_t N>
void printarray(std::array<int, N>& onehot) {
    printf("[ ");
    for (int i = 0; i < N; i++) {
        printf("%d ", onehot[i]);
    }
    printf("] \n");
}

int main() {
    Molecule mol;
    Array<char> descr;
    Array<char> symbol;
    loadMolecule("C1ON(CC)C1", mol);
    std::array<int, ONEHOT_SIZE> onehot;

    // loadMolecule(CAFFEINE, mol);
    for (int v : mol.vertices()) {
        int a = mol.getAtomNumber(v);
        mol.getAtomSymbol(v, symbol);
        mol.getAtomDescription(v, descr);

        // std::cout << atom_num << std::string(descr.ptr(), descr.size() - 1) << std::endl;
        onehot_encode(a, onehot);
        // printf("v: %d %s %d\n", v, Element::toString(a), a);

        bool in_ring = mol.vertexInRing(v);
        bool is_arom = mol.getAtomAromaticity(v) == ATOM_AROMATIC;

        printf("Atom    : %s\n", Element::toString(a));
        printf("In ring : %d\n", in_ring);
        printf("Aromatic: %d\n", is_arom);

        printarray(onehot);


    }

    std::cout << smiles(mol) << std::endl;
    return 0;
}
