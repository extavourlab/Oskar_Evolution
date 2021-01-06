'''
Authors: Savandara Besse and Leo Blondel

Creation: 07-05-2017
Modificaion: 09-02-2019

Description:
- Required inputs:
    - PDB structures
        - 5A4A.pdb (OSK)
        - 5NT7.pdb (LOTUS-VASA)
        - 2DB3.pdb (VASAR-NA)
    - CSV tables
        - ./scores.csv (all conservation scores saved in a CSV file generated with scripts in generate_Scores)
        - ./mapping.csv (mapping table between oskar alignment positions and oskar structure positions)

- Guidelines:
    - Run this script in Pymol
    - Command line template: SHOW, STRUCTURENAME, SCORE, AREA (option)
        - STRUCTURENAME types:Hydro
            - LOTUSVASA
            - OSK
        - SCORE types : Raw score - LOTUS monomers - LOTUS dimers - Holometabola - Hemimetabola
            - (Conservation) Valdar, Valdarmon, Valdardim, Valdarholo, Valdarhemi
            - (Conservation) JSD, JSDmon, JSDdim, JSDholo, JSDhemi
            - (Electrostaticity) Elec, Elecmon, Elecdim, Elecholo, Elechemi
            - (Hydrophobicity) Hydro, Hydromon, Hydrodim, Hydroholo, Hydrohemi
            - (Valdar Holo / Hemi) Logratio
            - (RNA-binding) RNABindR
        - AREA types :
            - For OSK --> OSK, core, rna-binding
            - For LOTUSVASA --> LOTUS, dimerization-interface, vasa-interface
    - LOTUSRNA function to generate RNA / LOTUSVASA alignment

'''


from pymol import cmd, stored
import pandas as pd
import numpy as np
import random


class helpers:
    def __init__(self, mapping, scores):
        self.mapping = pd.read_csv(mapping)
        self.scores = pd.read_csv(scores)

    def map(self, residue, domain):
        domain = self.mapping[self.mapping['domain']==domain]
        resi = domain[domain['structure'] == residue]
        return resi['alignment'].values[0]

    def score(self, score_type, position):
        score = self.scores[self.scores['position']==position]
        s = score[score_type].values[0]
        return s

    def get_score_by_resi(self, resi, domain, score_type):
        posi = self.map(resi, domain)
        score = self.score(score_type, posi)
        return score

def setup_pymol():
    cmd.set('scene_buttons', 1)
    cmd.set('matrix_mode', 1)
    cmd.set('movie_panel', 1)

    cmd.mset('1 x500')   #Number of frames
    cmd.set('ray_trace_frames', 1)
    cmd.viewport(800, 800)  #Resolution (px)

def init_representation(name):
    global lastpdbID
    if 'OSK' == name:
        pdbID = '5A4A'
    if 'LOTUSVASA' == name:
        pdbID = '5NT7'
    if 'VASA' == name :
        pdbID = '2DB3'
    saveview = True
    if lastpdbID != pdbID:
        cmd.fetch(pdbID)  # Load the PDB file
        cmd.hide('everything', 'all')  # Hide everything
        cmd.show('cartoon', 'all')  # Protein in cartoon
        if pdbID == '5A4A':
            cmd.hide('everything', 'resi 399-401')
            cmd.hide('everything', 'resi 604-606')
        lastpdbID = pdbID
        saveview = False

    return pdbID, saveview

def score(name,pdbID,function,score_type):
    scores = []
    other_Osk = []
    other_Lotus = []
    structure = {
        'OSK':'402,606',
        'LOTUSVASA':'146,238',
        'LOTUS':'148,224'
    }

    conserved_Osk = {
        427:'OSK',
        428:'OSK',
        429:'OSK',
        430:'OSK',
        431:'OSK',
        432:'OSK',
        424:'OSK',
        435:'OSK',
        438:'OSK',
        441:'OSK',
        452:'OSK',
        457:'OSK',
        561:'OSK',
        566:'OSK',
        593:'OSK',
        595:'OSK',
        436:'OSK',
        442:'OSK',
        576:'OSK',
        429:'OSK',
        458:'OSK',
        486:'OSK',
        487:'OSK',
        571:'OSK',
        586:'OSK'
    }

    conserved_Lotus = {
        155:'LOTUS',
        158:'LOTUS',
        159:'LOTUS',
        161:'LOTUS',
        162:'LOTUS',
        165:'LOTUS',
        167:'LOTUS',
        184:'LOTUS',
        224:'LOTUS',
        227:'LOTUS',
        228:'LOTUS',
        231:'LOTUS',
        235:'LOTUS',
        197:'LOTUS',
        200:'LOTUS',
        210:'LOTUS',
        215:'LOTUS'
    }

    colors = {
        'Valdar':'cyan_white_magenta',
        'JSD':'cyan_white_magenta',
        'Logratio':'cyan_white_magenta',
        'Hydro':'cyan_white_green',
        'Elec':'red_white_blue',
        'RNABindR':'yellow_white_red'
    }

    for key in colors.keys():
        if key in score_type:
            color = colors[key]

    Range = structure[name].split(',')
    for i in range(int(Range[0]),int(Range[1])):
        score = H.get_score_by_resi(i, name, score_type)
        cmd.alter("resi %s"%i, "b=%s"%score)
        scores.append(score)

        if 'OSK' in function:
            if i not in conserved_Osk.keys():
                other_Osk.append(str(i))
        if 'LOTUS' in function:
            if i not in conserved_Lotus.keys():
                other_Lotus.append(str(i))

    if 'OSK' in function:
        important_Osk = '+'.join([str(key) for key in conserved_Osk.keys()])
        cmd.select(function,'resi {}'.format(important_Osk))
    if 'LOTUS' in function:
        important_Lotus = '+'.join([str(key) for key in conserved_Lotus.keys()])
        cmd.select(function,'resi {}'.format(important_Lotus))

    if 'Elec' in score_type:
        cmd.spectrum("b",color, pdbID, -1, 1)
    elif 'Hydro' in score_type:
        cmd.spectrum("b",color, pdbID, -1, 1)
    elif 'Logratio' in score_type:
        cmd.spectrum("b",color, pdbID, 0, 2)
    else:
        cmd.spectrum("b",color, pdbID, 0, 1)
        if 'LOTUSVASA' in name and 'RNABindR' in score_type:
            cmd.select('LOTUS-1','/5NT7//A')
            cmd.show("mesh", 'LOTUS-1')
            cmd.select('LOTUS-2','/5NT7//C')
            cmd.show("mesh", 'LOTUS-2')
        if 'OSK' in name and 'RNABindR' in score_type:
            cmd.show("surface", pdbID)

    if 'OSK' in function:
        cmd.show('sticks',function)
        cmd.select('other_Osk','resi {}'.format('+'.join(other_Osk)))
        cmd.color('white','other_Osk')
    if 'LOTUS' in function:
        cmd.show('sticks',function)
        cmd.select('other_Lotus','resi {}'.format('+'.join(other_Lotus)))
        cmd.color('white','other_Osk')

    cmd.deselect()
    cmd.zoom(pdbID)

    if 'LOTUSVASA' in name :
        cmd.color("paleyellow", '/%s//D'%pdbID)
        cmd.color("paleyellow", '/%s//B'%pdbID)
        cmd.select('lotus-interface','resi 460+461+462+463+464+465+466+467+468+469+504+508+527')
        cmd.show('sticks','lotus-interface')
        cmd.select('rna-interface','resi 525+528')
        cmd.show('sticks','rna-interface')
        cmd.deselect()

    if function.strip() != '':
        focus_on(pdbID,function)

def focus_on(pdbID,function):
    toHighlight = {
        155:'vasa-interface',
        158:'vasa-interface',
        159:'vasa-interface',
        161:'vasa-interface',
        162:'vasa-interface',
        165:'vasa-interface',
        167:'vasa-interface',
        184:'vasa-interface',
        224:'vasa-interface',
        227:'vasa-interface',
        228:'vasa-interface',
        231:'vasa-interface',
        235:'vasa-interface',
        197:'dimerization-interface',
        200:'dimerization-interface',
        210:'dimerization-interface',
        215:'dimerization-interface',
        427:'core',
        428:'core',
        429:'core',
        430:'core',
        431:'core',
        436:'rna-binding',
        442:'rna-binding',
        576:'rna-binding'
    }

    if function != 'None':
        List = []
        for position in toHighlight.keys():
            if function == toHighlight[position]:
                List.append(str(position))
        nbRange = '+'.join(List)
        cmd.select(str(function),'resi {}'.format(nbRange))
        cmd.show('sticks',function)
        adjust_view(function)


def adjust_view(function):
    setview = {
        'OSK': '0.857618511, -0.507815480, -0.081356436, -0.421081513, -0.602513075, -0.677992046, 0.295276403, 0.615710795, -0.730552316, 0.000000000, 0.000000000, -120.490249634,    37.067352295,   20.989597321,   45.526634216,    82.197814941,  158.782684326,  -20.000000000',
        'core': '-0.616111696, -0.761785686, 0.200225979, -0.149057508, 0.362367362, 0.920033813, -0.773429453, 0.536997616, -0.336807847, -0.000373862, -0.000103349, -36.626907349, 37.511211395, 19.299640656, 40.756092072, -1.664040923, 74.920852661, -20.000000000',
        'rna-binding': '0.875713170, 0.415968120, 0.245161340, 0.445433080, -0.500040710, -0.742669821, -0.186332628, 0.759559333, -0.623170495, 0.000379004, -0.000125569,  -77.244255066,    34.662864685,   22.570341110,   43.849826813,    38.962089539,  115.546958923,  -20.000000000',
        'LOTUS': '0.347249299,    0.439891487,    0.828193963,   -0.236469552,   -0.813530445,    0.531257749,    0.907467067,   -0.380327016,   -0.178472668,    0.000203600,    0.000689294, -103.510437012,  -26.995281219,   34.309635162,   76.606254578,   44.202167511,  162.890518188,  -20.000000000',
        'dimerization-interface': '0.447004229,   -0.040959846,   -0.893581450,    -0.770385146,    0.490019232,   -0.407842278,     0.454581380,    0.870725691,    0.187487051,     0.000049177,   -0.003517030,  -43.662452698,   -28.647928238,   33.160560608,   80.513961792,   -16.673383713,  102.014999390,  -20.000024796',
        'vasa-interface': '0.145184740,   -0.984607041,    0.097281128, -0.565786779,   -0.163278759,   -0.808210313, 0.811651230,    0.062305421,   -0.580797255, -0.000039861, 0.000486811,  -85.625335693, -26.916032791,   48.033569336,   65.723213196, 26.272987366,  144.961410522,  -20.000000000'}

    closeup = setview[function]
    cmd.set_view(closeup)
    cmd.deselect()

def show(name,score_type,function=''):
    global lastpdbID
    pdbID, saveview = init_representation(name)
    if saveview:
        lastview = cmd.get_view()
    score(name,pdbID,function,score_type)
    if saveview:
        cmd.set_view(lastview)

def LOTUSRNA():
    cmd.fetch('5NT7')
    cmd.fetch('2DB3')
    cmd.hide('everything', '2DB3')
    cmd.show('cartoon', '2DB3')

    cmd.select('VASA','/2DB3//D')
    cmd.color("lightpink", 'VASA')
    cmd.select('RNA','/2DB3//H/')
    cmd.show('sticks','RNA')

    cmd.cealign('5NT7','2DB3', object='alignment')

    toRemove = {
        '2DB3':['A','B','C','E','F','G']
    }

    for i in range(len(toRemove['2DB3'])) :
        cmd.remove('/{}//{}/'.format('2DB3',toRemove['2DB3'][i]))

    show('LOTUSVASA','RNABindR','lotus-rna')

def save(name, res=1000):
    cmd.bg_color('white')
    cmd.set('ray_trace_mode', 0)
    cmd.set('antialias', 2)
    cmd.ray(res, res)
    cmd.png(name,res, res, res, 1)
    print('Picture saved!')

def save_all(name, filename, res=1000):
    score_types = ['Valdar', 'JSD', 'Elec', 'Hydro']
    for score_type in score_types:
        types = ['', 'holo', 'hemi', 'mon', 'dim']
        for type in types:
            show(name, score_type+type)
            save(filename + '_' + score_type + '_' + type + '.png', res)
    for score_type in ['Logratio', 'RNABindR']:
        show(name, score_type)
        save(filename + '_' + score_type + '.png', res)

mapping_table = '../Data/03_Oskar_scores_generation/CSV/mapping.csv'
score_table = '../Data/03_Oskar_scores_generation/CSV/scores.csv'

H = helpers(mapping_table, score_table)

lastpdbID = ""
cmd.extend('SHOW',show)
cmd.extend('LOTUSRNA',LOTUSRNA)
cmd.extend('SAVE',save)
cmd.extend('SAVEALL',save_all)
