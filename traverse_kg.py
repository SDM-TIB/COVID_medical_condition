import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON


def get_ddi(endpoint, comorb_drug, cov_drug):
    sparql = SPARQLWrapper(endpoint)
    list_drug = comorb_drug + cov_drug
    input_db_uri = ','.join(['<http://covid-19.tib.eu/vocab/'+db+'>' for db in list_drug])

    query = """select distinct ?drugLabel1 ?drugLabel2 ?effectLabel
                WHERE {
                ?DrugDrugInteraction a <http://covid-19.tib.eu/vocab/DrugDrugInteraction>.
                #?DrugDrugInteraction <http://covid-19.tib.eu/vocab/hasImpact> ?Impact.
                ?DrugDrugInteraction <http://covid-19.tib.eu/vocab/Effect> ?effect.
                ?DrugDrugInteraction <http://covid-19.tib.eu/vocab/precipitant> ?Drug1.
                ?DrugDrugInteraction <http://covid-19.tib.eu/vocab/object> ?Drug2.
                ?DrugDrugInteraction <http://covid-19.tib.eu/vocab/hasDrug1CUI> ?Drug1CUI.
                ?Drug1CUI <http://covid-19.tib.eu/vocab/drugLabel> ?drugLabel1.

                ?DrugDrugInteraction <http://covid-19.tib.eu/vocab/hasDrug2CUI> ?Drug2CUI.        
                ?Drug2CUI <http://covid-19.tib.eu/vocab/drugLabel> ?drugLabel2.
                ?effect <http://covid-19.tib.eu/vocab/drugLabel> ?effectLabel.

                FILTER (?drugLabel1 != ?drugLabel2)
                FILTER (?Drug1 in (""" + input_db_uri + """))

                FILTER (?Drug2 in (""" + input_db_uri + """))
    }"""

    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    dd = {'drugLabel1':[], 'drugLabel2':[], 'effectLabel':[]}
    for r in results['results']['bindings']:
        dd['drugLabel1'].append(r['drugLabel1']['value'])
        dd['drugLabel2'].append(r['drugLabel2']['value'])
        dd['effectLabel'].append(r['effectLabel']['value'])

    set_DDIs = pd.DataFrame(dd)
    #set_DDIs['Effect'] = set_DDIs['Effect'].str.replace('http://covid-19.tib.eu/Effect/', '')
    
    return set_DDIs