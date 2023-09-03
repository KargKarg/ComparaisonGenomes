import random
import distance
import numpy as np


def consensus(sequence1: str, sequence2: str) -> str:
    """
    Fonction de calcul de la séquence consensus entre 2 séquences.
    Un nucléotide sera préféré face à un gap.
    L'aléatoire choisira l'un des deux nucléotides s'ils sont différents.

    Argument:
        - sequence1 (str): La première séquence.
        - sequence2 (str): La deuxième séquence.

    Return:
         - cons (str): La séquence consensus.
    """
    cons = ""
    for i in range(len(sequence1)):
        if sequence1[i] == sequence2[i]:
            cons += sequence1[i]
        elif sequence1[i] == '-':
            cons += sequence2[i]
        elif sequence2[i] == '-':
            cons += sequence1[i]
        else:
            cons += random.choice([sequence1[i], sequence2[i]])

    return cons


def neighbor_joining(data: dict) -> str:
    """
    Fonction de calcul de la phylogénie entre 3 coregenomes avec la méthode NJ (neighbor joining).
    Calcul en premier temps la matrice des distances.

    Puis la difference net:
        - U(A) = Somme(D(A, G)) pour tout G appartenant aux Génomes étudiés.

    Puis calcul la matrice net:
        - D*(A, B) = D(A, B) - U(A) - U(B)

    Puis la distance entre les feuilles et le noeud K:
        - D(A, K) = (D(A, B) + U(A) - U(B)) / 2

    Les résultats sont disponibles dans le répertoire Neighbor Joining.
    Ils incluent la matrice des distances, la matrice net et enfin la phylogénie associée

    Argument:
        - data (dict): Dictionnaire sous la forme Nom:CGaligné.

    Return:
         - phylo (str): La phylogénie au format NEWICK.
    """
    annuaire = {indice: genome for indice, genome in enumerate(data.keys())}
    matrice_dist = np.full((len(data), len(data)), 0, float)
    matrice_net = np.full((len(data), len(data)), 0, float)
    U = [0 for _ in range(len(annuaire))]
    minimum = 1000000
    mi, mj = 0, 0

    for i in range(matrice_dist.shape[0]):
        for j in range(i+1, matrice_dist.shape[1]):
            matrice_dist[i, j] = distance.jaccard(100, data[annuaire[i]], data[annuaire[j]])
            matrice_dist[j, i] = matrice_dist[i, j].copy()
        U[i] = sum(matrice_dist[i])/(len(annuaire)-2)

    for i in range(matrice_net.shape[0]):
        for j in range(i+1, matrice_net.shape[1]):
            matrice_net[i, j] = matrice_dist[i, j] - U[i] - U[j]
            matrice_net[j, i] = matrice_net[i, j].copy()
            if matrice_net[i, j] < minimum:
                mi, mj = i, j
                minimum = matrice_net[i, j]

    cons = consensus(data[annuaire[mi]], data[annuaire[mj]])
    phylo = f"({annuaire[mi]}:{round((matrice_dist[mi, mj]+U[mi]-U[mj])/2, 4)}, {annuaire[mj]}:{round((matrice_dist[mi, mj]+U[mj]-U[mi])/2, 4)})"

    del data[annuaire[mi]], data[annuaire[mj]]

    val = round(distance.jaccard(100, cons, list(data.values())[0]), 4)

    phylo = f"({list(data.keys())[0]}:{val/2}, {phylo}:{val/2})"

    return phylo


def upgma(data: dict) -> str:
    """
    Fonction de calcul de la phylogénie entre 3 coregenomes avec la méthode UPGMA.
    Calcul en premier temps la matrice des distances.
    Puis calcul la distance du dernier:
        - D(C, (A,B)) = D(C, A) * 0.5 + D(C, B) * 0.5

    Les résultats sont disponibles dans le répertoire UPGMA.
    Ils incluent la matrice des distances et la phylogénie associée

    Argument:
        - data (dict): Dictionnaire sous la forme Nom:CGaligné.

    Return:
         - phylo (str): La phylogénie au format NEWICK.
    """
    phylo = ""
    annuaire = {indice: genome for indice, genome in enumerate(data.keys())}
    annuaire_inverse = {genome: indice for indice, genome in enumerate(data.keys())}
    matrice_dist = np.full((len(data), len(data)), 0, float)
    minimum = 1000000
    mi, mj = 0, 1

    for i in range(matrice_dist.shape[0]):
        for j in range(i+1, matrice_dist.shape[1]):
            matrice_dist[i, j] = distance.jaccard(100, data[annuaire[i]], data[annuaire[j]])
            matrice_dist[j, i] = matrice_dist[i, j].copy()
            if matrice_dist[i, j] < minimum:
                mi, mj = i, j
                minimum = matrice_dist[i, j]

    phylo += f"({annuaire[mi]}:{round(matrice_dist[mi, mj] / 2, 4)}, {annuaire[mj]}:{round(matrice_dist[mi, mj] / 2, 4)})"
    data = list(data.keys())
    data.remove(annuaire[mi])
    data.remove(annuaire[mj])
    val = round((matrice_dist[annuaire_inverse[data[0]], mi] + matrice_dist[annuaire_inverse[data[0]], mj])/2, 4)
    phylo = f"({data[0]}:{val}, {phylo}:{val})"

    return phylo
