import json
from utils import eq_affiliation, eq_name



class Publications:

    INSTITUTIONS_TEXAS = [
        'University of Texas Rio Grande Valley',
        'West Texas A&M University',
        'University of North Texas Health Science Center',
        'Texas A&M University-Commerce',
        'Texas A&M University',
        'Southern Methodist University',
        'Texas A&M International University',
        'Sam Houston State University',
        'University of St. Thomas',
        'University of North Texas',
        'University of Texas at San Antonio',
        'University of Texas Southwestern Medical Center',
        'University of Texas at Arlington',
        'Texas State University',
        'Texas Christian University',
        'Baylor University',
        'University of Texas Medical Branch at Galveston',
        'University of Texas Medical Branch, Galveston',
        'Texas A&M University-Kingsville',
        'Prairie View A&M University',
        'University of Texas at Dallas',
        'Tarleton State University',
        'East Texas A&M University',
        'Texas Tech University Health Sciences Center',
        'Texas Tech University',
        'Trinity University',
        "St. Mary's University",
        'University of Texas at Austin',
        'Abilene Christian University',
        'Texas A&M University Health Science Center',
        'University of Houston-Clear Lake',
        'Rice University',
        'Baylor College of Medicine',
        'University of Texas Permian Basin',
        'University of Texas at El Paso',
        'University of Houston',
        'Midwestern State University',
        'University of Texas M. D. Anderson Cancer Center',
        'Lamar University',
        "Texas Woman's University",
        'University of Texas Health Science Center at Houston'
    ]
    
    @staticmethod
    def filter_publications():
        publications = Publications.load_from_json("data/publications.json") or []
        db_researchers = Publications.load_from_json("data/researchers.json") or []
        pubs_in_file = set(p.get('doi') for p in Publications.load_from_json("data/filtered_publications.json"))
        filtered_pubs = []
        print(len(publications), len(db_researchers), len(pubs_in_file))
        for pub in publications:
            if pub.get('doi') in pubs_in_file:
                continue
            pub_crossref = pub.get('crossref')
            if pub_crossref.get('access'):
                pub_checked = Publications.flag_authors(pub_crossref, db_researchers)
                if pub_checked:
                    pub['crossref'] = pub_checked
                    filtered_pubs.append(pub)
        Publications.add_publication_to_json(filtered_pubs, "data/filtered_publications.json")

    @staticmethod
    def generate_summaries_for_publications():
        from llm_tools import generate_summary
        publications = Publications.load_from_json("data/filtered_publications.json") or []
        pubs_display = Publications.load_from_json("data/publications_display.json") or []
        pubs_with_summaries_dois = set(p.get('doi') for p in pubs_display)
        for pub in publications:
            if pub.get('doi') in pubs_with_summaries_dois:
                continue
            title = pub.get('title', '')
            abstract = pub.get('crossref', {}).get('abstract')
            if not abstract:
                pubs_display.append(pub)
                continue
            summary = generate_summary(title, abstract)
            if summary:
                pub['summary'] = summary
                pubs_display.append(pub)
        Publications.add_publication_to_json(pubs_display, "data/publications_display.json")

    @staticmethod
    def load_from_json(filename):
        try:
            with open(filename, "r") as f:
                return json.load(f)
        except:
            return None
    
    @staticmethod
    def add_publication_to_json(data, filename):
        existing_data = Publications.load_from_json(filename) or []
        for pub in data:
            if not any(p.get('doi') == pub.get('doi') for p in existing_data):
                existing_data.append(pub)
        with open(filename, "w") as f:
            json.dump(existing_data, f, indent=2)

    @staticmethod
    def flag_authors(publication, db_researchers):
        """Check if any author of the publication is affiliated with a Texas institution."""
        include_pub = False
        for author in publication.get('authors', []):
            first_name = author.get('first_name', '')
            last_name = author.get('last_name', '')
            affiliations = author.get('institution', [])
            orcid = author.get('orcid')
            for affiliation in affiliations:
                for texas_inst in Publications.INSTITUTIONS_TEXAS:
                    if eq_affiliation(texas_inst, affiliation) or 'texas' in affiliation.lower():
                        author['texas_affiliation'] = True
                        include_pub = True
                        break
            for researcher in db_researchers:
                if orcid and researcher.get('orcid_id') and orcid == researcher.get('orcid_id'):
                    author['db_researcher_id'] = researcher.get('researcher_id') 
                    include_pub = True 
                    break
                elif last_name == researcher.get('last_name') and researcher.get('first_name') == first_name:
                    if len(affiliations) > 0:
                        for affiliation in affiliations:
                            if eq_affiliation(researcher.get('institutions'), affiliation):
                                author['db_researcher_id'] = researcher.get('researcher_id')
                                include_pub = True
                                break
                        break
                    else: 
                        author['db_researcher_id'] = researcher.get('researcher_id')
                        break
        if include_pub:
            return publication
        return None

def fix_publications():
    from journals import Journal
    publications = Publications.load_from_json("data/publications_display.json") or []
    # publications_filter = Publications.load_from_json("data/filtered_publications.json") or []
    for pub in publications:
        print(pub['crossref'].get('published_date'), pub.get('published_date'))
        p = Journal.get_publication_crossref(pub['doi'])
        p_date = p['published_date']
        if len(p_date) < 3:
            rss_date = pub.get('published_date', '')
            p_date = [int(x) for x in rss_date.split('-')]
            print(f"Using RSS date for {pub['doi']}: {p_date}")
        print(pub['doi'], p_date)
        pub['crossref']['published_date'] = p_date
        
    save_path = "data/publications_display_fixed.json"
    with open(save_path, "w") as f:
        json.dump(publications, f, indent=2)

# fix_publications()

# Publications.filter_publications()
# Publications.generate_summaries_for_publications()

# print(len(Publications.INSTITUTIONS_TEXAS))
