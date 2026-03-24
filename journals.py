import requests
import time
import json
from Bio import Entrez
from rss_parsers import RSSParserFactory


Entrez.email = "vito.rozman@revelo.bi"
Entrez.api_key = '833bee09582b616923d32401ca51cede1508'


class Publication:
    """Fetch and parse publication metadata from Crossref and PubMed."""

    @staticmethod
    def _get_abstract(doi, crossref_data):
        """Extract abstract from Crossref or PubMed."""
        pub = Publication.from_pubmed_doi(doi)
        if pub.get('access') is None:
            return pub.get('abstract', '').strip().replace('\n', ' ')
        abstract = crossref_data.get('abstract', '')
        return abstract.strip().replace('\n', ' ') if abstract else None

    @staticmethod
    def from_crossref(doi):
        """Get publication metadata from Crossref API by DOI."""
        response = requests.get(f"https://api.crossref.org/works/{doi}", timeout=30)
        if response.status_code != 200:
            return {'access': False, 'doi': doi}
        data = response.json().get('message', {})
        return {
            'access': True,
            'doi': doi,
            'title': data.get('title', [''])[0],
            'published_date': data.get('published', {}).get('date-parts', [[]])[0],
            'authors': [{
                'first_name': a.get("given", ""),
                'last_name': a.get("family", ""),
                'institution': [aff.get('name') for aff in a.get("affiliation", []) if 'name' in aff],
                'orcid': a.get("ORCID", "").replace("https://orcid.org/", "") if a.get("ORCID") else None
            } for a in data.get("author", [])],
            'abstract': Publication._get_abstract(doi, data)
        }

    @staticmethod
    def from_pubmed_doi(doi):
        """Get publication from PubMed by DOI."""
        for attempt in range(3):
            try:
                handle = Entrez.esearch(db="pubmed", term=f"{doi}[AID]", retmax=1)
                results = Entrez.read(handle)
                handle.close()
                if results["IdList"]:
                    return Publication.from_pubmed_pmids(results["IdList"])[0]
                return {'access': False, 'doi': doi}
            except Exception:
                if attempt < 2:
                    time.sleep(2 ** (attempt + 1))
        return {'access': False, 'doi': doi}

    @staticmethod
    def from_pubmed_pmids(pmids):
        """Get publications from PubMed by PMIDs."""
        months = {'Jan':'01','Feb':'02','Mar':'03','Apr':'04','May':'05','Jun':'06',
                  'Jul':'07','Aug':'08','Sep':'09','Oct':'10','Nov':'11','Dec':'12'}
        for attempt in range(3):
            try:
                handle = Entrez.efetch(db="pubmed", id=",".join(pmids), rettype="gb", retmode="xml")
                records = Entrez.read(handle)
                handle.close()
                break
            except Exception:
                if attempt < 2:
                    time.sleep(2 ** (attempt + 1))
                else:
                    return []
        
        publications = []
        for record in records["PubmedArticle"]:
            article = record["MedlineCitation"]["Article"]
            doi = next((str(id_) for id_ in record["PubmedData"]["ArticleIdList"] 
                       if id_.attributes["IdType"] == "doi"), None)
            date_data = (article.get("ArticleDate", [{}])[0] if article.get("ArticleDate") 
                        else article["Journal"]["JournalIssue"]["PubDate"])
            yyyy, mm, dd = date_data.get("Year"), date_data.get("Month", "01"), date_data.get("Day", "01")
            date = f"{yyyy}-{months.get(mm, str(mm).zfill(2))}-{str(dd).zfill(2)}" if yyyy else None
            authors = [{
                'first_name': a.get("ForeName"),
                'last_name': a.get("LastName"),
                'institution': [aff.get("Affiliation") for aff in a.get("AffiliationInfo", []) 
                              if "Affiliation" in aff]
            } for a in article.get("AuthorList", []) if a.get("LastName") and a.get("ForeName")]
            
            publications.append({
                "doi": doi,
                "title": article.get("ArticleTitle", ""),
                "abstract": str(article.get("Abstract", {}).get("AbstractText", [""])[0]),
                "date": date,
                "authors": authors
            })
        return publications


class Journal:
    """Fetch publications from journal RSS feeds and enrich with metadata."""
    FEEDS_FILE = "data/journal_rss_feeds.json"

    @staticmethod
    def _resolve_feeds_file(feeds_file=None):
        return feeds_file or Journal.FEEDS_FILE

    @staticmethod
    def get_journal_feeds(feeds_file=None):
        path = Journal._resolve_feeds_file(feeds_file)
        records = Journal.load_from_json(path) or []
        return {
            item.get("name"): item.get("rss_link")
            for item in records
            if item.get("name") and item.get("rss_link")
        }

    @staticmethod
    def get_new_publications(journal, feed_url):
        parser = RSSParserFactory.get_parser(journal, feed_url)
        return parser.parse(feed_url)

    @staticmethod
    def load_from_json(filename):
        try:
            with open(filename, "r") as f:
                return json.load(f)
        except Exception:
            return None
    
    @staticmethod
    def add_publication_to_json(data, filename):
        existing_data = Journal.load_from_json(filename) or []
        for pub in data:
            if not any(p.get('doi') == pub.get('doi') for p in existing_data):
                existing_data.append(pub)
        with open(filename, "w") as f:
            json.dump(existing_data, f, indent=2)




    @staticmethod
    def fetch_publications_from_feeds(file_path='data/publications.json', feeds_file=None, verbose=True):
        dois_in_file = set(p.get('doi') for p in Journal.load_from_json(file_path) or [])
        all_publications = []
        feeds = Journal.get_journal_feeds(feeds_file)
        for journal, feed_url in feeds.items():
            publications = []
            if verbose:
                print(f"Fetching publications from {journal}...")
            pubs_info = Journal.get_new_publications(journal, feed_url)
            for pub in pubs_info:
                if not pub['doi'] or pub['doi'] in dois_in_file:
                    continue
                item = {
                    'journal_name': journal,
                    'doi': pub['doi'],
                    'title': pub['title'],
                    'link': pub['link'],
                    'published_date': pub['published_date'],
                    'crossref': Publication.from_crossref(pub['doi']),
                }
                publications.append(item)
                dois_in_file.add(pub['doi'])
                time.sleep(0.3)  # Sleep to avoid hitting API rate limits
            if verbose:
                print(f"  Found {len(pubs_info)} publications, {len(publications)} new.")
            all_publications.extend(publications)
            Journal.add_publication_to_json(publications, file_path)
        return all_publications


class JournalTesting:
    """Small helper class for manual parsing/metadata checks."""

    @staticmethod
    def _get_feeds_for_testing(journal_name=None, feeds_file=None):
        feeds = Journal.get_journal_feeds(feeds_file)
        if not journal_name:
            return feeds

        exact_name = journal_name
        if exact_name not in feeds:
            exact_name = {name.lower(): name for name in feeds}.get(journal_name.lower())

        if not exact_name:
            print(f"Journal '{journal_name}' not found in feeds file.")
            return None

        return {exact_name: feeds[exact_name]}

    @staticmethod
    def test_first_publication_parsing(journal_name=None, feeds_file=None):
        """Print first parsed publication for one journal or all feeds."""
        feeds_to_test = JournalTesting._get_feeds_for_testing(journal_name, feeds_file)
        if feeds_to_test is None:
            return {}

        results = {}
        for journal, feed_url in feeds_to_test.items():
            
            try:
                parsed_items = Journal.get_new_publications(journal, feed_url)
                if not parsed_items:
                    print(f"\n=== {journal} ===")
                    print("No entries parsed.")
                    results[journal] = None
                    continue

                first_publication = parsed_items[0]
                # print only if doi is null:
                if not first_publication.get('doi'):
                    print(f"\n=== {journal} ===")
                    print(json.dumps(first_publication, indent=2, ensure_ascii=False))
                results[journal] = first_publication
            except Exception as exc:
                print(f"Error while parsing: {exc}")
                results[journal] = None

        return results

    @staticmethod
    def run_manual_checks():
        """Run default manual checks."""
        JournalTesting.test_first_publication_parsing()

    @staticmethod
    def test_first_publication_details(journal_name=None, feeds_file=None):
        """Fetch Crossref details for the first parsed publication with a valid DOI."""

        feeds_to_test = JournalTesting._get_feeds_for_testing(journal_name, feeds_file)
        if feeds_to_test is None:
            return {}

        results = {}
        for journal, feed_url in feeds_to_test.items():
            
            try:
                parsed_items = Journal.get_new_publications(journal, feed_url)
                first_with_doi = next((item for item in parsed_items if item.get('doi')), None)
                if not first_with_doi:
                    print(f"\n=== {journal} ===")
                    print("No publication with valid DOI found.")
                    results[journal] = None
                    continue

                doi = first_with_doi['doi']
                details = Publication.from_crossref(doi)

                date = "-".join([str(i) for i in details.get('published_date')]) if details.get('published_date') else "N/A"
                abstract_info = "Has abstract" if details.get('abstract') else "No abstract"
                if details.get('access') is False:
                    continue
                elif abstract_info == "No abstract":
                    print(f"\n=== {journal} ===")
                    print(f"DOI {details.get('doi') == doi}: {doi}")
                    print(f"Date: {date}, {abstract_info}")
                # print(json.dumps(details, indent=2, ensure_ascii=False))

                results[journal] = {
                    'doi': doi,
                    'details': details,
                }
            except Exception as exc:
                print(f"Error while processing {journal}: {exc}")
                results[journal] = None

        return results
    

########################################################################################
# Testing

if __name__ == "__main__":
    Journal.fetch_publications_from_feeds(verbose=True)
    # JournalTesting.run_manual_checks()