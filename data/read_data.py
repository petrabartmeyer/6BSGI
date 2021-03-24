from pandas import read_csv
from os.path import isdir, join


def load_data_as_json(dataset, output='output.json'):

    p = load_data(dataset)

    from json import dump

    with open(output, 'w') as fp:

        dump(p, fp)


def load_data(dataset):
    """This function reads all the files of a dataset (name of a
    directory) and returns a dictionary containing all the
    information. If the dataset does not exist, returns an empty
    dictionary.

    """

    if not isdir(dataset):

        return {}

    usinas = read_csv(join(dataset, 'info.csv'), index_col=0, skipinitialspace=True)

    # Create problem
    problem = {i: {'nUG': usinas.loc[i, 'NUG']} for i in usinas.index}
    # Create units
    for i in problem.keys():
        problem[i]['UG'] = list({} for j in range(problem[i]['nUG']))

    # Add general information
    cmon = read_csv(join(dataset, 'cota_montante.csv'), index_col=0, sep="\s*,\s*", engine='python')
    for i in problem.keys():
        problem[i]['cotaMontante'] = cmon.loc[i].to_list()

    cjus = read_csv(join(dataset, 'cota_jusante.csv'), index_col=0, sep="\s*,\s*", engine='python')
    for i in problem.keys():
        problem[i]['cotaJusante'] = cjus.loc[i].to_list()

    lim = read_csv(join(dataset, 'limites.csv'), index_col=0, sep="\s*,\s*", engine='python')
    for u in lim.index:
        problem[u].update(lim.loc[u])

    # Add information to each unit
    add_dict_to_unity(problem, join(dataset, 'limites_potencia.csv'))
    add_dict_to_unity(problem, join(dataset, 'rampa.csv'))

    add_list_to_unity(problem, 'vazaoMax', join(dataset, 'vazao_turbinada_maxima.csv'))
    add_list_to_unity(problem, 'vazaoMin', join(dataset, 'vazao_turbinada_minima.csv'))
    add_list_to_unity(problem, 'perdaGerador', join(dataset, 'perda_gerador.csv'))
    add_list_to_unity(problem, 'perdaHidraulica', join(dataset, 'perda_hidraulica.csv'))
    add_list_to_unity(problem, 'perdaMecTurbina', join(dataset, 'perda_mecanica_turbina.csv'))
    add_list_to_unity(problem, 'rendimentoHidraulico', join(dataset, 'rendimento_hidraulico.csv'))
    
    return problem


def add_dict_to_unity(problem, fname):
    """
    This function updates the dictionary of units.
    """

    pot = read_csv(fname, index_col=[0, 1], sep="\s*,\s*", engine='python')
    for u, i in pot.index:
        ugs = problem[u]['UG']
        if i is -1:
            for j in range(problem[u]['nUG']):
                ugs[j].update(pot.loc[(u, i)])
        else:
            ugs[i].update(pot.loc[(u, i)])


def add_list_to_unity(problem, key, fname):
    """
    This function adds a list (usually parameters of a polynomial) to the unit.
    """

    df = read_csv(fname, index_col=[0, 1], sep="\s*,\s*", engine='python')
    for u, i in df.index:
        ugs = problem[u]['UG']
        if i is -1:
            for j in range(problem[u]['nUG']):
                ugs[j][key] = df.loc[(u, i)].to_list()
        else:
            ugs[i][key] = df.loc[(u, i)].to_list()
