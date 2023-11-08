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
    with open('CoreGenome/core_genome.txt', 'r') as fil:
        for ligne in fil:
            donnees_cg[ligne.replace('\n', '')[1:]] = fil.__next__().replace('\n', '')

    return donnees_cg


def core_genome_align() -> dict:
    """
    Fonction qui prend les CG alignés sous format FASTA et qui crée un dictionnaire avec clé:valeur

    clé: Nom du génome.
    valeur: CG aligné.

    Argument:
        - None

    Return:
        - donnees_cg_align (dict): structure contenant les CG alignés.
    """
    donnees_cg_align = {}
    with open('CoreGenome/core_genome_alignement.txt', 'r') as fil:
        for ligne in fil:
            donnees_cg_align[ligne.replace('\n', '')[1:]] = fil.__next__().replace('\n', '')

    return donnees_cg_align


def caracteres() -> dict:
    """
    Fonction qui prend les CG alignés en dictionnaire avec clé:valeur et qui séléctionne les caractères informatifs
    dans le but de réaliser une phylogénie.

    Ces caractères ont:
        - au moins 1 état différent
        - ne présente pas de gap
        - au moins 2 états identiques

    clé: Nom du génome.
    valeur: CG aligné.

    Argument:
        - None

    Return:
        - carac (dict): structure contenant les caractères.
    """
    cg = core_genome_align()
    carac = {gid: "" for gid in cg.keys()}
    sequences = list(cg.values())
    for i in range(len(sequences[0])):
        etats = set()
        for seq in sequences:
            etats = etats.union(seq[i])
        if '-' not in etats and len(etats) != 1:
            for gid, seq in cg.items():
                carac[gid] += f"{seq[i]}"

    with open('CoreGenome/caracteres.txt', 'w') as filout:
        texte = ""
        for gid, car in carac.items():
            texte += f">{gid}\n{car}\n"
        filout.write(texte)

    return carac
