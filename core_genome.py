import alignement_pair
import alignement_triple
import extraction
import threading


def clustering(g1: list, g2: list, g3: list, n1: str, n2: str, n3: str, tppfid: int) -> None:
    """
    Fonction qui crée des clusters de 3 protéines si le %id est supérieur à 94%.
    Utilise l'alignement global de Needleman et Wunsch pour déduire un alignement global de score maximum.
    Si Alignement(seq1, seq2) > 94%id et Alignement(seq1, seq3) > 94%id et Alignement(seq2, seq3) > 94%id:
        Alors les 3 protéines sont considérés dans le CoreGenome.

    Les résultats sont disponibles dans le dossier Resultats/coregenome_clusters.txt, les cds de gauches appartiennent
    au Génome 1, celles du milieu au Génome 2 et celles de gauche au Génome 3.

    Un fichier texte comprenant le Core Génome du Génome spécifique est aussi créée.

    Argument:
        - g1 (list): contient les CDS du Génome 1.
        - g2 (list): contient les CDS du Génome 2.
        - g3 (list): contient les CDS du Génome 3.
        - n1 (str): nom du Génome 1.
        - n2 (str): nom du Génome 2.
        - n3 (str): nom du Génome 3.
        - tpfid (int): ID du thread petit fils.

    Return:
        - None.
    """
    cpt1 = 0
    for seq1 in g1:
        cpt1 += 1
        cpt2 = 0
        for seq2 in g2:
            cpt2 += 1
            cpt3 = 0
            trouver = False

            if alignement_pair.pourcentage_identite(seq1, seq2) > 94:

                for seq3 in g3:
                    cpt3 += 1

                    with open(f"Rapports/threads_{tppfid}.txt", 'a') as fil:
                        fil.write(f'{cpt1}/{len(g1)}   {cpt2}/{len(g2)}   {cpt3}/{len(g3)}\n')

                    if alignement_pair.pourcentage_identite(seq1, seq3) > 94 and alignement_pair.pourcentage_identite(seq2, seq3) > 94:
                        with open(f"Resultats/core_genome_{n1}.txt", 'a') as fil:
                            fil.write(f"{seq1}\n")
                        with open(f"Resultats/core_genome_{n2}.txt", 'a') as fil:
                            fil.write(f"{seq2}\n")
                        with open(f"Resultats/core_genome_{n3}.txt", 'a') as fil:
                            fil.write(f"{seq3}\n")
                        with open("Resultats/core_genome_clusters.txt", 'a') as fil:
                            fil.write(f"{seq1};{seq2};{seq3}\n")
                        trouver = True
                        break

                if trouver:
                    break

            else:
                with open(f"Rapports/threads_{tppfid}.txt", 'a') as fil:
                    fil.write(f'{cpt1}/{len(g1)}   {cpt2}/{len(g2)}   {cpt3}/{len(g3)}\n')


def petitpetitfils(g1: list, g2: list, g3: list, n1: str, n2: str, n3: str, tpfid: int) -> None:
    """
    Fonction qui crée des threads petit petit fils pour paralléliser le clustering.
    Chaque threads prendra une partie des CDS du Génome 1 pour les comparer avec une partie du Génome 2 et une partie
    du Génome 3.

    Argument:
        - g1 (list): contient les CDS du Génome 1.
        - g2 (list): contient les CDS du Génome 2.
        - g3 (list): contient les CDS du Génome 3.
        - n1 (str): nom du Génome 1.
        - n2 (str): nom du Génome 2.
        - n3 (str): nom du Génome 3.
        - tfid (int): ID du thread petit fils.

    Return:
        - None.
    """
    nb_petitpetitfils = 4

    sous_tab = len(g3) // nb_petitpetitfils
    restant = len(g3) % nb_petitpetitfils

    t1 = threading.Thread()

    cpt = 0

    for i in range(0, len(g3), sous_tab):
        cpt += 1
        if i == len(g3) - sous_tab - restant:
            t1 = threading.Thread(target=clustering, args=(g1, g2, g3[i:], n1, n2, n3, tpfid+cpt))
        else:
            t1 = threading.Thread(target=clustering, args=(g1, g2, g3[i:i+sous_tab], n1, n2, n3, tpfid+cpt))
        t1.start()
    t1.join()


def petitfils(g1: list, g2: list, g3: list, n1: str, n2: str, n3: str, tfid: int) -> None:
    """
    Fonction qui crée des threads petit fils pour paralléliser le clustering.
    Chaque threads prendra une partie des CDS du Génome 1 pour les comparer avec une partie du Génome 2 et l'ensemble
    du Génome 3.

    Argument:
        - g1 (list): contient les CDS du Génome 1.
        - g2 (list): contient les CDS du Génome 2.
        - g3 (list): contient les CDS du Génome 3.
        - n1 (str): nom du Génome 1.
        - n2 (str): nom du Génome 2.
        - n3 (str): nom du Génome 3.
        - tfid (int): ID du thread fils.

    Return:
        - None.
    """
    nb_petitfils = 4

    sous_tab = len(g2) // nb_petitfils
    restant = len(g2) % nb_petitfils

    t1 = threading.Thread()

    cpt = 0

    for i in range(0, len(g2), sous_tab):
        cpt += 1
        if i == len(g2) - sous_tab - restant:
            t1 = threading.Thread(target=petitpetitfils, args=(g1, g2[i:], g3, n1, n2, n3, tfid+cpt*10))
        else:
            t1 = threading.Thread(target=petitpetitfils, args=(g1, g2[i:i+sous_tab], g3, n1, n2, n3, tfid+cpt*10))
        t1.start()
    t1.join()


def fils() -> None:
    """
    Fonction qui crée des threads pour paralléliser le clustering.
    Chaque threads prendra une partie des CDS du Génome 1 pour les comparer avec l'ensemble des Génome 2 et Génome 3.

    Argument:
        - None.

    Return:
        - None.
    """
    nb_fils = 4

    cds = extraction.cds()

    n1, n2, n3 = cds.keys()
    g1, g2, g3 = cds.values()

    sous_tab = len(g1) // nb_fils
    restant = len(g1) % nb_fils

    t1 = threading.Thread()

    cpt = 0

    for i in range(0, len(g1), sous_tab):
        cpt += 1
        if i == len(g1) - sous_tab - restant:
            t1 = threading.Thread(target=petitfils, args=(g1[i:], g2, g3, n1, n2, n3, cpt*100))
        else:
            t1 = threading.Thread(target=petitfils, args=(g1[i:i + sous_tab], g2, g3, n1, n2, n3, cpt*100))
        t1.start()
    t1.join()


def mise_en_forme(n1: str, n2: str, n3: str) -> None:
    """
    Fonction qui concatène les séquences codantes pour former une séquence continue.
    Bien que l'ordre des CDS soit aléatoire, chaque homologue est positionné au même niveau que les autres.

    Argument:
        - n1 (str): Nom du Génome 1.
        - n2 (str): Nom du Génome 2.
        - n3 (str): Nom du Génome 3.

    Return:
        - None.
    """
    with open(f"Resultats/core_genome_{n1}.txt", 'r') as fil:
        seq1 = ""
        for ligne in fil:
            seq1 += ligne.replace('\n', '')

    with open(f"Resultats/core_genome_{n2}.txt", 'r') as fil:
        seq2 = ""
        for ligne in fil:
            seq2 += ligne.replace('\n', '')

    with open(f"Resultats/core_genome_{n3}.txt", 'r') as fil:
        seq3 = ""
        for ligne in fil:
            seq3 += ligne.replace('\n', '')

    with open("Resultats/core_genome.txt", 'w') as fil:
        fil.write(f">{n1}\n{seq1}\n>{n2}\n{seq2}\n>{n3}\n{seq3}")


def alignement() -> None:
    """
    Fonction qui aligne les 3 CG et les sauvegarde dans Resultats au format FASTA.

    Argument:
        - None.

    Return:
        - None.
    """
    n1, n2, n3 = extraction.core_genome().keys()
    cg1, cg2, cg3 = extraction.core_genome().values()

    align1, align2, align3 = alignement_triple.alignement_global(cg1, cg2, cg3)

    with open('Resultats/core_genome_alignement.txt', 'w') as fil:
        fil.write(f">{n1}\n{align2}\n>{n2}\n{align2}\n>{n3}\n{align3}")
