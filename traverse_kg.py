import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON


def store_pharmacokinetic_ddi(effect):
    if effect in ['excretion_rate', 'excretory_function', 'excretion', 'excretion rate']:
        effect = 'excretion'
    elif effect in ['process_of_absorption', 'process of absorption', 'absorption']:
        effect = 'absorption'
    elif effect in ['serum_concentration', 'serum_concentration_of', 'serum_level', 'serum_globulin_level',
                    'metabolite', 'active_metabolites', 'serum concentration']:
        effect = 'serum_concentration'
    elif effect in ['metabolism']:
        effect = 'metabolism'
    else:
        if effect == 'adverse effects':
            effect = 'adverse_effects'
        elif effect == 'corrected_prolonged_qt_interval_by_ecg_finding':
            effect = 'prolonged_qt'
        elif effect == 'corrected prolonged qt interval by ecg finding':
            effect = 'prolonged_qt'
        return effect, True
    return effect, False


def get_ddi(endpoint, comorb_drug, cov_drug):
    sparql = SPARQLWrapper(endpoint)
    list_drug = comorb_drug + cov_drug
    input_db_uri = ','.join(['<http://research.tib.eu/covid-19/entity/'+db+'>' for db in list_drug])
    
    query = """
    prefix covid-19: <http://research.tib.eu/covid-19/vocab/>
    select distinct ?precipitant_label ?object_label ?effect_label ?impact ?precipitant_db ?object_db
                        where {

                        ?DrugDrugInteraction a covid-19:DrugDrugInteraction.        
                        ?DrugDrugInteraction covid-19:precipitant_hasDrugBankID ?precipitant_db.
                        ?precipitant_db covid-19:hasCUIAnnotation ?precipitant_cui.
                        ?precipitant_cui covid-19:annLabel ?precipitant_label.
                        ?DrugDrugInteraction covid-19:object_hasDrugBankID ?object_db.
                        ?object_db covid-19:hasCUIAnnotation ?object_cui.
                        ?object_cui covid-19:annLabel ?object_label.
                        ?DrugDrugInteraction covid-19:effect ?effect.
                        ?effect covid-19:annLabel ?effect_label .
                        ?DrugDrugInteraction covid-19:impact ?impact.

                        FILTER (?precipitant_label != ?object_label)
                        FILTER (?precipitant_db in (""" + input_db_uri + """))
                        FILTER (?object_db in (""" + input_db_uri + """))
                       
                        }"""

    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    dd = {'precipitant_label':[], 'object_label':[], 'effect_label':[], 'impact':[], 'precipitant_db':[], 'object_db':[]}
    for r in results['results']['bindings']:
        dd['precipitant_label'].append(r['precipitant_label']['value'].lower())
        dd['object_label'].append(r['object_label']['value'].lower())
        effect = r['effect_label']['value'].lower()
        effect, pharmadynamic = store_pharmacokinetic_ddi(effect)
        dd['effect_label'].append(effect.lower())
        dd['impact'].append(r['impact']['value'])
        dd['precipitant_db'].append(r['precipitant_db']['value'])
        dd['object_db'].append(r['object_db']['value'])

    set_DDIs = pd.DataFrame(dd)
    
    cov_drug_uri = ['http://research.tib.eu/covid-19/entity/'+db+'' for db in cov_drug]
    a = set(set_DDIs.precipitant_db.unique())
    a.update(list(set_DDIs.object_db.unique()))
    b = a.intersection(set(cov_drug_uri))
    cov_ddi = [w.replace('http://research.tib.eu/covid-19/entity/', '') for w in b]
    
    # set_DDIs['effect_label'] = set_DDIs['effect_label'].str.replace('corrected_prolonged_qt_interval_by_ecg_finding', 'prolonged_qt')
    set_DDIs.drop(set_DDIs[set_DDIs.precipitant_label=='injection,_azithromycin,_500_mg_administered'].index, inplace=True)
    set_DDIs.drop(set_DDIs[set_DDIs.object_label=='injection,_azithromycin,_500_mg_administered'].index, inplace=True)


    set_DDIs['effect_impact'] = set_DDIs[['effect_label', 'impact']].apply(lambda x: '_'.join(x), axis=1)
    set_DDIs = set_DDIs[['precipitant_label', 'object_label', 'effect_impact']]
    set_DDIs.drop_duplicates(inplace=True)
    return set_DDIs, cov_ddi

def get_drug_label(endpoint, cov_drug):
    sparql = SPARQLWrapper(endpoint)
    input_db_uri = ['<http://research.tib.eu/covid-19/entity/'+db+'> ' for db in cov_drug]
    
    drugLabel = []
    for d in input_db_uri:
        query = """
        prefix covid-19: <http://research.tib.eu/covid-19/vocab/>
        select distinct ?drugLabel
                     where {
                        """+ d +""" covid-19:hasCUIAnnotation ?drug_cui.
                        ?drug_cui covid-19:annLabel ?drugLabel.
                    }"""
        sparql.setQuery(query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()

        for r in results['results']['bindings']:
            drugLabel.append(r['drugLabel']['value'].lower())
    return drugLabel