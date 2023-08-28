import math
import numpy as np


def hamming(sequence1: str, sequence2: str) -> int:
    """
    Fonction de calcul de la distance Hamming entre deux séquences.

    Argument:
        - sequence1 (str): La première séquence.
        - sequence2 (str): La deuxième séquence.

    Return:
         - (int): La distance Hamming entre les deux séquences.
    """
    dst = 0
    for i in range(len(sequence1)):
        if sequence1[i] != sequence2[i]:
            dst += 1

    return dst


def hamming_revisite(sequence1: str, sequence2: str) -> int:
    """
    Fonction de calcul de la distance Hamming entre deux séquences.

    En se basant sur la probabilité accrue des transitions,
    la distance augmente plus ou moins en fonction de la mutation associée.

    Argument:
        - sequence1 (str): La première séquence.
        - sequence2 (str): La deuxième séquence.

    Return:
         - (int): La distance Hamming entre les deux séquences.
    """
    dst = 0
    purine, pyrimidine = {'A', 'G'}, {'C', 'T'}

    for i in range(len(sequence1)):
        if sequence1[i] != sequence2[i]:
            if sequence1[i] in purine and sequence2[i] in purine:
                dst += 1
            elif sequence1[i] in pyrimidine and sequence2[i] in pyrimidine:
                dst += 1
            else:
                dst += 2

    return dst


def kmer(k: int, sequence: str) -> set:
    """
    Fonction de calcul de l'ensemble des kmers d'une séquence.

    Argument:
        - k (int): Taille des mers considérés.

    Return:
         - (set): L'ensemble des kmers de la séquence.
    """
    mer = set()
    for i in range(0, len(sequence)-k+1):
        mer.add(sequence[i:i+k])

    return mer


def jaccard(k: int, sequence1: str, sequence2: str) -> float:
    """
    Fonction de calcul de la distance Jaccard entre deux séquences.
    Calcul en premier temps les ensembles de kmer propre à chaques séquences.
    Puis calcul la distance Jaccard grâce à la formule:
        - (1-|AUB|/|A&B|)

    Argument:
        - k (int): Taille des mers considérés.
        - sequence1 (str): La première séquence.
        - sequence2 (str): La deuxième séquence.

    Return:
         - (float): La distance Jaccard entre les deux séquences.
    """
    mer1 = kmer(k, sequence1)
    mer2 = kmer(k, sequence2)

    return 1 - (len(mer1.intersection(mer2))/len(mer1.union(mer2)))


def jukes_cantor(sequence1: str, sequence2: str) -> float:
    """
    Fonction de calcul de la distance Jukes-Cantor entre deux séquences.
    Calcul en premier temps la proportion de nt différent entre les deux séqunences gràce à la distance Hamming.
    Puis calcul la distance JC grâce à la formule:
        - (-3/4)*(ln(1-4/3*p)
        p: La proportion de nt différent

    Argument:
        - sequence1 (str): La première séquence.
        - sequence2 (str): La deuxième séquence.

    Return:
         - (float): La distance Jukes-Cantor entre les deux séquences.
    """
    p = hamming(sequence1, sequence2) / len(sequence1)

    return (-3/4)*(math.log(1-(4/3)*p))


def levenshtein(sequence1: str, sequence2: str) -> int:
    """
    Fonction de calcul de la distance Levenshtein entre deux séquences.

    Argument:
        - sequence1 (str): La première séquence.
        - sequence2 (str): La deuxième séquence.

    Return:
         - (int): La distance Levenshtein entre les deux séquences.
    """
    matrice = np.full((len(sequence1)+1, len(sequence2)+1), 0)

    for i in range(1, matrice.shape[0]):
        matrice[i, 0] = matrice[i-1, 0] + 1

    for i in range(1, matrice.shape[1]):
        matrice[0, i] = matrice[0, i-1] + 1

    for i in range(1, matrice.shape[0]):
        for j in range(1, matrice.shape[1]):
            if sequence1[i-1] == sequence2[j-1]:
                matrice[i, j] = matrice[i-1, j-1]
            else:
                matrice[i, j] = 1 + min(matrice[i-1, j], matrice[i, j-1], matrice[i-1, j-1])

    return matrice[matrice.shape[0]-1, matrice.shape[1]-1]
