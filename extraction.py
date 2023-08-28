import os


def cds() -> dict:
    """
    Fonction qui prend les CDS sous format FASTA et qui crée un dictionnaire avec clé:valeur

    clé: Nom du génome.
    valeur: liste contenant l'ensemble des CDS.

    Argument:
        - None

    Return:
        - donnees_cds (dict): structure contenant les CDS.
    """
    donnee_cds = {}
    for strain in os.listdir('Donnees'):
        donnee_cds[strain] = []
        with open(f"Donnees/{strain}/cds.fna", 'r') as filin:
            coding = ""
            for ligne in filin:
                if ligne[0] == ">":
                    donnee_cds[strain].append(coding)
                    coding = ""
                else:
                    coding += ligne.replace('\n', '')
            donnee_cds[strain].append(coding)
        donnee_cds[strain] = donnee_cds[strain][1:]

    return donnee_cds


def core_genome() -> dict:
    """
    Fonction qui prend les CG sous format FASTA et qui crée un dictionnaire avec clé:valeur

    clé: Nom du génome.
    valeur: CG.

    Argument:
        - None

    Return:
        - donnees_cg (dict): structure contenant les CG.
    """
    donnees_cg = {}
    with open('Resultats/core_genome.txt', 'r') as fil:
        for ligne in fil:
            donnees_cg[ligne.replace('\n', '')[1:]] = fil.__next__().replace('\n', '')

    return donnees_cg
