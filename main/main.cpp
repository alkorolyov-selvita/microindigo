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

#include "molecule/smiles_loader.h"


using namespace indigo;

void loadMolecule(const char* buf, Molecule& molecule)
{
    BufferScanner scanner(buf);
    // MoleculeAutoLoader loader(scanner);
    SmilesLoader loader(scanner);
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

    map[ELEM_C] = 0;
    map[ELEM_N] = 1;
    map[ELEM_O] = 2;
    map[ELEM_P] = 3;
    map[ELEM_S] = 4;
    map[ELEM_F] = 5;
    map[ELEM_Cl] = 6;
    map[ELEM_Br] = 7;
    map[ELEM_I] = 8;

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
    AromaticityOptions options(AromaticityOptions::GENERIC);

    bool in_ring, is_aromatic, is_hbd, is_hba;
    int hetero_nei_count, heavy_nei_count;

    std::array<int, ONEHOT_SIZE> onehot;

    // loadMolecule("C1ON(CC)C1", mol);
    loadMolecule(CAFFEINE, mol);
    mol.aromatize(options);

    std::cout << "Atom Ring HANei HetNei Arom Donor OneHot" << std::endl;

    for (int v : mol.vertices()) {
        int a = mol.getAtomNumber(v);

        // mol.getAtomSymbol(v, symbol);
        // mol.getAtomDescription(v, descr);
        // std::cout << atom_num << std::string(descr.ptr(), descr.size() - 1) << std::endl;
        // printf("v: %d %s %d\n", v, Element::toString(a), a);


        onehot_encode(a, onehot);
        in_ring = mol.vertexInRing(v);
        is_aromatic = mol.getAtomAromaticity(v) == ATOM_AROMATIC;

        hetero_nei_count = 0;
        heavy_nei_count = 0;

        const Vertex& vertex = mol.getVertex(v);
        int i;
        int idx;
        int nei;
//        printf("vertex: %d\n", v);
        for (i = vertex.neiBegin(); i != vertex.neiEnd(); i = vertex.neiNext(i)) {
            idx = vertex.neiVertex(i);
//            printf("  idx, nei: %d %d\n", i, nei);
            nei = mol.getAtomNumber(idx);
            if (nei != a) hetero_nei_count ++;
            if (nei != ELEM_H) heavy_nei_count ++;
        }

        // HBD SMARTS
        // "N&!H0&v3 N&!H0&+1&v4
        // O&H1&+0 S&H1&+0
        // n&H1&+0
        // ver 2.0.1
        int h_count = mol.getAtomTotalH(v);
        int charge = mol.getAtomCharge(v);
        int valence = mol.getAtomValence(v);

        bool nitrogen_aliphatic_donor =
            a == ELEM_N &&
            h_count != 0 &&
            (
                charge == 0 && valence == 3 ||
                charge == 1 && valence == 4
            );
        bool oxygen_donor = a == ELEM_O && h_count == 1 && charge == 0;
        bool sulphur_donor = a == ELEM_S && h_count == 1 && charge == 0;
        bool nitrogen_aromatic_donor = a == ELEM_N && is_aromatic && h_count == 1 && charge == 0;

        bool is_donor = nitrogen_aliphatic_donor || oxygen_donor || sulphur_donor || nitrogen_aromatic_donor;


        // HBA SMARTS
        // $([O,S;H1;v2]-[!$(*=[O,N,P,S])])
        // $([O,S;H0;v2]) $([O,S;-])
        // $([N;v3;!$(N-*=!@[O,N,P,S])])
        // $([nH0,o,s;+0])
        // ver 2.0.1
        bool hba_conditions[4] = {false};

        // $([O,S;H1;v2]-[!$(*=[O,N,P,S])])
        hba_conditions[0] = (a == ELEM_O || a == ELEM_S) && h_count == 1 && valence == 2;
        // check nei
        bool nei_condition = false;
        for (i = vertex.neiBegin(); i != vertex.neiEnd(); i = vertex.neiNext(i)) {
            idx = vertex.neiVertex(i);
            nei = mol.getAtomNumber(idx);

        }




        printf("[%s]: %d %d %d %d %d ", Element::toString(a), in_ring, heavy_nei_count, hetero_nei_count, is_aromatic, is_donor);
        printarray(onehot);

    }

    std::cout << smiles(mol) << std::endl;
    return 0;
}
