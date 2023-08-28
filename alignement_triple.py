import numpy as np


def lecture_triple(matrice: np.array, sequence1: str, sequence2: str, sequence3: str, sid: int, ssub: int, sgap: int) -> tuple:
    """
    Fonction qui permet la lecture de la  matrice 3D avec l'algorithme de Needleman et Wunsch.
    Renvoie les trois séquences alignées.

    Arguments:
        - matrice (np.array) Matrice remplie.
        - sequence1 (str): Première séquence à aligner.
        - sequence2 (str): Deuxième séquence à aligner.
        - sequence3 (str): Troisième séquence à aligner.
        - sid (int): Score d'une identité.
        - ssub (int): Score d'une substitution.
        - sgap (int): Score d'un gap.

    Return:
        - (tuple): Les trois séquences alignées.
    """
    align1, align2, align3 = "", "", ""
    i, j, k = len(sequence1), len(sequence2), len(sequence3)

    while i > 0 or j > 0 or k > 0:

        voisin1 = matrice[i - 1, j - 1, k - 1]
        voisin1 += sid if sequence1[i - 1] == sequence2[j - 1] == sequence3[k - 1] else ssub
        voisin2 = matrice[i - 1, j - 1, k] + sgap
        voisin2 += sid if sequence1[i - 1] == sequence2[j - 1] else ssub
        voisin3 = matrice[i - 1, j, k - 1] + sgap
        voisin3 += sid if sequence1[i - 1] == sequence3[k - 1] else ssub
        voisin4 = matrice[i, j - 1, k - 1] + sgap
        voisin4 += sid if sequence2[j - 1] == sequence3[k - 1] else ssub
        voisin5 = matrice[i - 1, j, k] + 2 * sgap
        voisin6 = matrice[i, j - 1, k] + 2 * sgap
        voisin7 = matrice[i, j, k - 1] + 2 * sgap

        if i > 0 and j > 0 and k > 0:
            suivant = max(voisin1, voisin2, voisin3, voisin4, voisin5, voisin6, voisin7)

            if suivant == voisin1:
                align1, align2, align3 = sequence1[i - 1] + align1, sequence2[j - 1] + align2, sequence3[k - 1] + align3
                i, j, k = i - 1, j - 1, k - 1
            elif suivant == voisin2:
                align1, align2, align3 = sequence1[i - 1] + align1, sequence2[j - 1] + align2, '-' + align3
                i, j = i - 1, j - 1
            elif suivant == voisin3:
                align1, align2, align3 = sequence1[i - 1] + align1, '-' + align2, sequence3[k - 1] + align3
                i, k = i - 1, k - 1
            elif suivant == voisin4:
                align1, align2, align3 = '-' + align1, sequence2[j - 1] + align2, sequence3[k - 1] + align3
                j, k = j - 1, k - 1
            elif suivant == voisin5:
                align1, align2, align3 = sequence1[i - 1] + align1, '-' + align2, '-' + align3
                i -= 1
            elif suivant == voisin6:
                align1, align2, align3 = '-' + align1, sequence2[j - 1] + align2, '-' + align3
                j -= 1
            else:
                align1, align2, align3 = '-' + align1, '-' + align2, sequence3[k - 1] + align3
                k -= 1

        elif i > 0 and k > 0:
            suivant = max(voisin3, voisin5, voisin7)

            if suivant == voisin3:
                align1, align2, align3 = sequence1[i - 1] + align1, '-' + align2, sequence3[k - 1] + align3
                i, k = i - 1, k - 1
            elif suivant == voisin5:
                align1, align2, align3 = sequence1[i - 1] + align1, '-' + align2, '-' + align3
                i -= 1
            else:
                align1, align2, align3 = '-' + align1, '-' + align2, sequence3[k - 1] + align3
                k -= 1

        elif i > 0 and j > 0:
            suivant = max(voisin2, voisin5, voisin6)

            if suivant == voisin2:
                align1, align2, align3 = sequence1[i - 1] + align1, sequence2[j - 1] + align2, '-' + align3
                i, j = i - 1, j - 1
            elif suivant == voisin5:
                align1, align2, align3 = sequence1[i - 1] + align1, '-' + align2, '-' + align3
                i -= 1
            elif suivant == voisin6:
                align1, align2, align3 = '-' + align1, sequence2[j - 1] + align2, '-' + align3
                j -= 1

        elif j > 0 and k > 0:
            suivant = max(voisin4, voisin6, voisin7)

            if suivant == voisin4:
                align1, align2, align3 = '-' + align1, sequence2[j - 1] + align2, sequence3[k - 1] + align3
                j, k = j - 1, k - 1
            elif suivant == voisin6:
                align1, align2, align3 = '-' + align1, sequence2[j - 1] + align2, '-' + align3
                j -= 1
            else:
                align1, align2, align3 = '-' + align1, '-' + align2, sequence3[k - 1] + align3
                k -= 1

        elif i > 0:
            align1, align2, align3 = sequence1[i - 1] + align1, '-' + align2, '-' + align3
            i -= 1

        elif j > 0:
            align1, align2, align3 = '-' + align1, sequence2[j - 1] + align2, '-' + align3
            j -= 1

        elif k > 0:
            align1, align2, align3 = '-' + align1, '-' + align2, sequence3[k - 1] + align3
            k -= 1

    return align1, align2, align3


def alignement_global(sequence1: str, sequence2: str, sequence3: str) -> tuple:
    """
    Fonction qui crée et initialise une matrice 3D avec l'algorithme de Needleman et Wunsch.
    Réalise ainsi un alignement global de 3 séquences simultanément.
    Renvoie les trois séquences alignées.
    Les paramètres sont fixes:
        - identité: +2
        - substitution: -1
        - gap: -2

    Arguments:
        - sequence1 (str): Première séquence à aligner.
        - sequence2 (str): Deuxième séquence à aligner.
        - sequence3 (str): Troisième séquence à aligner.

    Return:
        - (tuple): Les trois séquences alignées.
    """
    matrice = np.full((len(sequence1)+1, len(sequence2)+1, len(sequence3)+1), 0)

    sid, ssub, sgap = 2, -1, -2

    for i in range(1, len(sequence1)+1):
        matrice[i, 0, 0] = matrice[i-1, 0, 0] + sgap

    for j in range(1, len(sequence2)+1):
        matrice[0, j, 0] = matrice[0, j-1, 0] + sgap

    for k in range(1, len(sequence3)+1):
        matrice[0, 0, k] = matrice[0, 0, k-1] + sgap

    for i in range(1, len(sequence1)+1):
        for j in range(1, len(sequence2)+1):
            for k in range(1, len(sequence3)+1):
                voisin1 = matrice[i-1, j-1, k-1]
                voisin1 += sid if sequence1[i-1] == sequence2[j-1] == sequence3[k-1] else ssub
                voisin2 = matrice[i-1, j-1, k] + sgap
                voisin2 += sid if sequence1[i-1] == sequence2[j-1] else ssub
                voisin3 = matrice[i-1, j, k-1] + sgap
                voisin3 += sid if sequence1[i-1] == sequence3[k-1] else ssub
                voisin4 = matrice[i, j-1, k-1] + sgap
                voisin4 += sid if sequence2[j-1] == sequence3[k-1] else ssub
                voisin5 = matrice[i-1, j, k] + 2*sgap
                voisin6 = matrice[i, j-1, k] + 2*sgap
                voisin7 = matrice[i, j, k-1] + 2*sgap
                matrice[i, j, k] = max(voisin1, voisin2, voisin3, voisin4, voisin5, voisin6, voisin7)

    return lecture_triple(matrice, sequence1, sequence2, sequence3, sid, ssub, sgap)
