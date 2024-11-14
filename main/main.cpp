// #include "main.h"

#include "main.h"

#include <iostream>
#include <string>

#include "base_cpp/output.h"
#include <base_cpp/scanner.h>
#include <molecule/molecule.h>
#include "molecule/elements.h"
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
    int map[256] = {9}; // 9 for unknown atom

    // define known atoms
    map[ELEM_C] = 0;  // C
    map[ELEM_N] = 1;  // N
    map[ELEM_O] = 2;  // O
    map[ELEM_P] = 3;  // P
    map[ELEM_S] = 4;  // S
    map[ELEM_F] = 5;  // F
    map[ELEM_Cl] = 6; // Cl
    map[ELEM_Br] = 7; // Br
    map[ELEM_I] = 8;  // I

    int idx = map[atom_num];

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
    AromaticityOptions options;
    bool is_in_ring, is_aromatic, is_hbd, is_hba;
    int hetero_nei_count, heavy_nei_count;


    std::array<int, ONEHOT_SIZE> onehot;

    // loadMolecule("C1ON(CC)C1", mol);
    loadMolecule(CAFFEINE, mol);
    mol.aromatize(options);

    std::cout << "Atom Ring Arom Donor OneHot" << std::endl;

    for (int v : mol.vertices()) {
        int a = mol.getAtomNumber(v);

        // mol.getAtomSymbol(v, symbol);
        // mol.getAtomDescription(v, descr);
        // std::cout << atom_num << std::string(descr.ptr(), descr.size() - 1) << std::endl;
        // printf("v: %d %s %d\n", v, Element::toString(a), a);



        onehot_encode(a, onehot);
        bool in_ring = mol.vertexInRing(v);
        bool is_aromatic = mol.getAtomAromaticity(v) == ATOM_AROMATIC;


        // HBD SMARTS "[N&!H0&v3,N&!H0&+1&v4,O&H1&+0,S&H1&+0,n&H1&+0]" "2.0.1"
        int h_count = mol.getAtomTotalH(v);
        int charge = mol.getAtomCharge(v);
        int valence = mol.getAtomValence(v);

        bool is_nitrogen_aliphatic_hbd =
            a == ELEM_N &&
            h_count != 0 &&
            (
                charge == 0 && valence == 3 ||
                charge == 1 && valence == 4
            );
        bool is_oxygen_hbd = a == ELEM_O && h_count == 1 && charge == 0;
        bool is_sulphure_hbd = a == ELEM_S && h_count == 1 && charge == 0;
        bool is_nitrogen_aromatic_hbd = a == ELEM_N && is_aromatic && h_count == 1 && charge == 0;

        bool is_donor = is_nitrogen_aliphatic_hbd || is_oxygen_hbd || is_sulphure_hbd || is_nitrogen_aromatic_hbd;


        printf("[%s]: %d %d %d ", Element::toString(a), in_ring, is_aromatic, is_donor);
        printarray(onehot);


        // const Vertex& vertex = mol.getVertex(v);
        // int i;
        // is_in_ring = false;
        // is_aromatic = false;
        // for (i = vertex.neiBegin(); i != vertex.neiEnd(); i = vertex.neiNext(i)) {
        //     if (not is_in_ring & mol.getEdgeTopology(vertex.neiEdge(i)) == TOPOLOGY_RING) {
        //         is_in_ring = true;
        //     }
        //
        // }




    }

    std::cout << smiles(mol) << std::endl;
    return 0;
}
