import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON


def get_ddi(endpoint, comorb_drug, cov_drug):
    sparql = SPARQLWrapper(endpoint)
    list_drug = comorb_drug + cov_drug
    input_db_uri = ','.join(['<http://covid-19.tib.eu/Drug/'+db+'>' for db in list_drug])

    query = """select distinct ?precipitant_label ?object_label ?effect_label ?impact
                        where {

                        ?DrugDrugInteraction a <http://covid-19.tib.eu/vocab/DrugDrugInteraction>.        
                        ?DrugDrugInteraction <http://covid-19.tib.eu/vocab/precipitant_hasDrugBankID> ?precipitant_db.
                        ?precipitant_db <http://covid-19.tib.eu/vocab/hasCUIAnnotation> ?precipitant_cui.
                        ?precipitant_cui <http://covid-19.tib.eu/vocab/annlabel> ?precipitant_label.
                        ?DrugDrugInteraction <http://covid-19.tib.eu/vocab/object_hasDrugBankID> ?object_db.
                        ?object_db <http://covid-19.tib.eu/vocab/hasCUIAnnotation> ?object_cui.
                        ?object_cui <http://covid-19.tib.eu/vocab/annlabel> ?object_label.
                        ?DrugDrugInteraction <http://covid-19.tib.eu/vocab/effect> ?effect.
                        ?effect <http://covid-19.tib.eu/vocab/annlabel> ?effect_label .
                        ?DrugDrugInteraction <http://covid-19.tib.eu/vocab/impact> ?impact.

                        FILTER (?precipitant_label != ?object_label)
                        FILTER (?precipitant_db in (""" + input_db_uri + """))
                        FILTER (?object_db in (""" + input_db_uri + """))
                       
                        }"""

    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    dd = {'precipitant_label':[], 'object_label':[], 'effect_label':[], 'impact':[]}
    for r in results['results']['bindings']:
        dd['precipitant_label'].append(r['precipitant_label']['value'])
        dd['object_label'].append(r['object_label']['value'])
        dd['effect_label'].append(r['effect_label']['value'])
        dd['impact'].append(r['impact']['value'])

    set_DDIs = pd.DataFrame(dd)
    set_DDIs['effect_impact'] = set_DDIs[['effect_label', 'impact']].apply(lambda x: '_'.join(x), axis=1)
    set_DDIs = set_DDIs[['precipitant_label', 'object_label', 'effect_impact']]
    
    return set_DDIs

def get_drug_label(endpoint, cov_drug):
    sparql = SPARQLWrapper(endpoint)
    input_db_uri = ['<http://covid-19.tib.eu/Drug/'+db+'> ' for db in cov_drug]
    
    drugLabel = []
    for d in input_db_uri:
        query = """select distinct ?drugLabel
                     where {
                        """+ d +""" <http://covid-19.tib.eu/vocab/hasCUIAnnotation> ?drug_cui.
                        ?drug_cui <http://covid-19.tib.eu/vocab/annlabel> ?drugLabel.
                    }"""
        sparql.setQuery(query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()

        
        for r in results['results']['bindings']:
            drugLabel.append(r['drugLabel']['value'])
    
    return drugLabel
