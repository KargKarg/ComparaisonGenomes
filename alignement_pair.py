from Bio import pairwise2
from Bio.Seq import Seq


def pourcentage_identite(sequence1: str, sequence2: str) -> float:
    """
    Fonction qui crée et initialise une matrice 2D avec l'algorithme de Needleman et Wunsch.
    Renvoie le %id entre les 2 séquences.
    Les paramètres sont fixes:
        - identité: +1
        - substitution: 0
        - gap: 0

    Arguments:
        - sequence1 (str): Première séquence à aligner.
        - sequence2 (str): Deuxième séquence à aligner.

    Return:
        - pourc_id (float): Le pourcentage d'identité entre les 2 séquences.
    """
    alignments = pairwise2.align.globalxx(Seq(sequence1), Seq(sequence2))
    alignment = alignments[0]
    taille = max(len(sequence1), len(sequence2))
    id_match = alignment[-3]
    pourc_id = (id_match / taille) * 100
    return pourc_id


def sequences(sequence1: str, sequence2: str) -> tuple:
    """
    Fonction qui crée et initialise une matrice 2D avec l'algorithme de Needleman et Wunsch.
    Renvoie les deux séquences alignées.
    Les paramètres sont fixes:
        - identité: +1
        - substitution: 0
        - gap: 0

    Arguments:
        - sequence1 (str): Première séquence à aligner.
        - sequence2 (str): Deuxième séquence à aligner.

    Return:
        - (tuple): Les deux séquences alignées.
    """
    sequence1, sequence2 = Seq(sequence1), Seq(sequence2)
    return pairwise2.align.globalxx(sequence1, sequence2)[0][0], pairwise2.align.globalxx(sequence1, sequence2)[0][1]

