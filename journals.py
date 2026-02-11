import feedparser
import requests
from Bio import Entrez
import time
import json
from datetime import datetime

Entrez.email = "vito.rozman@revelo.bi"


class Journal:
    """Utility class for fetching scientific journal publications from various sources."""
    
    JOURNAL_FEEDS = {
        "Nature Chemistry": "https://www.nature.com/nchem.rss",
        "Nature Reviews Chemistry": "https://www.nature.com/natrevchem.rss",
        "Chemical Reviews": "http://pubs.acs.org/action/showFeed?jc=chreay&type=etoc&feed=rss",
        "Journal of the American Chemical Society": "http://pubs.acs.org/action/showFeed?jc=jacsat&type=etoc&feed=rss",
        "Angewandte Chemie Int. Ed.": "https://onlinelibrary.wiley.com/feed/15213773/most-recent",
        "Accounts of Chemical Research": "http://pubs.acs.org/action/showFeed?jc=achre4&type=etoc&feed=rss",
        "Advanced Materials": "https://onlinelibrary.wiley.com/feed/15214095/most-recent",
        "Nature Communications": "https://www.nature.com/ncomms.rss",
        "ACS Nano": "http://pubs.acs.org/action/showFeed?jc=ancac3&type=etoc&feed=rss",
        "Organic Letters": "http://pubs.acs.org/action/showFeed?jc=orlef7&type=etoc&feed=rss",
        "Inorganic Chemistry": "http://pubs.acs.org/action/showFeed?jc=inocaj&type=etoc&feed=rss",
        "Analytical Chemistry": "http://pubs.acs.org/action/showFeed?jc=ancham&type=etoc&feed=rss",
        "J. Phys. Chem. Letters": "http://pubs.acs.org/action/showFeed?jc=jpclcd&type=etoc&feed=rss",
    }

    @staticmethod
    def parse_date_from_feed(entry):
        """Parse publication date from RSS feed entry and return YYYY-MM-DD."""
        
        # Try ACS style
        prism_date = entry.get("prism_coverdate")
        if prism_date:
            return prism_date[:10]
        
        # Try Nature style
        updated_date = entry.get("updated")
        if updated_date and len(updated_date) >= 10:
            return updated_date[:10]
        
        # Try Wiley style
        published_str = entry.get("published")
        if published_str:
            try:
                return datetime.strptime(published_str, "%a, %d %b %Y %H:%M:%S %z").strftime("%Y-%m-%d")
            except (ValueError, TypeError):
                pass  # fail silently, will return None at the end
        
        return None

    @staticmethod
    def get_new_pubs_from_rss(feed_url, verify_ssl=True):
        """Extract DOIs from RSS feed URL."""
        feed = None
        try:
            response = requests.get(feed_url, timeout=10, headers={'User-Agent': 'Mozilla/5.0'}, verify=verify_ssl)
            response.raise_for_status()
            feed = feedparser.parse(response.content)
        except requests.exceptions.SSLError:
            if verify_ssl:
                try:
                    response = requests.get(feed_url, timeout=10, headers={'User-Agent': 'Mozilla/5.0'}, verify=False)
                    response.raise_for_status()
                    feed = feedparser.parse(response.content)
                except:
                    feed = None
        except:
            feed = None
        if not feed or 'entries' not in feed:
            return []
        pubs_info = []
        for entry in feed.get('entries', []):
            doi = entry.get("prism_doi")
            title = entry.get("title")
            link = entry.get("link")
            published_date = Journal.parse_date_from_feed(entry)
            pubs_info.append(
                {
                    "doi": doi,
                    "title": title,
                    "link": link,
                    "published_date": published_date
                })
        return pubs_info
    
    @staticmethod
    def get_dois_from_rss_rsc(feed_url, verify_ssl=True):
        """Extract DOIs from RSS feed URL."""
        import re
        
        def parse_doi_from_summary(summary):
            # Pattern matches DOI after "DOI:" or "<b>DOI</b>:" in various formats
            match = re.search(r'(?:DOI|<b>DOI</b>)\s*:\s*(10\.\d{4,}/[^\s,<]+)', summary, re.IGNORECASE)
            return match.group(1) if match else None
        
        feed = None
        try:
            response = requests.get(feed_url, timeout=10, headers={'User-Agent': 'Mozilla/5.0'}, verify=verify_ssl)
            response.raise_for_status()
            feed = feedparser.parse(response.content)
        except requests.exceptions.SSLError:
            if verify_ssl:
                try:
                    response = requests.get(feed_url, timeout=10, headers={'User-Agent': 'Mozilla/5.0'}, verify=False)
                    response.raise_for_status()
                    feed = feedparser.parse(response.content)
                except:
                    return []
        except:
            return []
        
        if not feed or 'entries' not in feed:
            return []
        
        dois = []
        for entry in feed.get('entries', []):
            doi = parse_doi_from_summary(entry.get("summary", ""))
            if doi:
                dois.append(doi)
        return dois
    
    @staticmethod
    def get_pubmed_ids(journal_query, retmax=100):
        """Get latest publication IDs from PubMed by journal query."""
        for attempt in range(3):
            try:
                handle = Entrez.esearch(db="pubmed", term=journal_query, sort="pub+date", retmax=retmax)
                results = Entrez.read(handle)
                handle.close()
                return results["IdList"]
            except Exception as e:
                if attempt < 2:
                    time.sleep(2 ** (attempt + 1))
        return []
    
    @staticmethod
    def get_publication_crossref(doi):
        """Get publication metadata from Crossref API by DOI."""
        response = requests.get(f"https://api.crossref.org/works/{doi}")
        if response.status_code != 200:
            return {'access': False, 'doi': doi}
        data = response.json().get('message', {})
        return {
            'access': True,
            'doi': doi,
            'title': data.get('title', [''])[0],
            'published_date': data.get('indexed', {}).get('date-parts', [[]])[0],
            'authors': [{
                'first_name': a.get("given", ""),
                'last_name': a.get("family", ""),
                'institution': [aff.get('name') for aff in a.get("affiliation", []) if 'name' in aff]
            } for a in data.get("author", [])]
        }
    
    @staticmethod
    def get_publication_pumbed_form_doi(doi):
        """Get publication metadata from PubMed by DOI."""
        for attempt in range(3):
            try:
                handle = Entrez.esearch(db="pubmed", term=f"{doi}[AID]", retmax=1)
                results = Entrez.read(handle)
                handle.close()
                if results["IdList"]:
                    return Journal.get_publications_pubmed(results["IdList"])[0]
                else:
                    return {'access': False, 'doi': doi}
            except Exception as e:
                if attempt < 2:
                    time.sleep(2 ** (attempt + 1))
        return {'access': False, 'doi': doi}

    @staticmethod
    def get_publications_pubmed(pmids):
        """Get publication metadata from PubMed by list of PMIDs."""
        months = {'Jan':'01','Feb':'02','Mar':'03','Apr':'04','May':'05','Jun':'06',
                  'Jul':'07','Aug':'08','Sep':'09','Oct':'10','Nov':'11','Dec':'12'}
        for attempt in range(3):
            try:
                handle = Entrez.efetch(db="pubmed", id=",".join(pmids), rettype="gb", retmode="xml")
                records = Entrez.read(handle)
                handle.close()
                break
            except:
                if attempt < 2:
                    time.sleep(2 ** (attempt + 1))
                else:
                    return []
        
        publications = []
        for record in records["PubmedArticle"]:
            article = record["MedlineCitation"]["Article"]
            doi = next((str(id_) for id_ in record["PubmedData"]["ArticleIdList"] if id_.attributes["IdType"] == "doi"), None)
            date_data = article.get("ArticleDate", [{}])[0] if article.get("ArticleDate") else article["Journal"]["JournalIssue"]["PubDate"]
            yyyy, mm, dd = date_data.get("Year"), date_data.get("Month", "01"), date_data.get("Day", "01")
            date = f"{yyyy}-{months.get(mm, str(mm).zfill(2))}-{str(dd).zfill(2)}" if yyyy else None
            authors = [{
                'first_name': a.get("ForeName"),
                'last_name': a.get("LastName"),
                'institution': [aff.get("Affiliation") for aff in a.get("AffiliationInfo", []) if "Affiliation" in aff]
            } for a in article.get("AuthorList", []) if a.get("LastName") and a.get("ForeName")]
            
            publications.append({
                "doi": doi,
                "title": article.get("ArticleTitle", ""),
                "abstract": str(article.get("Abstract", {}).get("AbstractText", [""])[0]),
                "date": date,
                "authors": authors
            })
        return publications

    @staticmethod
    def load_from_json(filename):
        try:
            with open(filename, "r") as f:
                return json.load(f)
        except:
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
    def fetch_publications_from_feeds(verbose=True):
        publications = []
        dois_in_file = set(p.get('doi') for p in Journal.load_from_json("data/publications.json") or [])
        for journal, feed_url in Journal.JOURNAL_FEEDS.items():
            if verbose:
                print(f"Fetching publications from {journal}...")
            pubs_info = Journal.get_new_pubs_from_rss(feed_url)
            for pub in pubs_info:
                if not pub['doi'] or pub['doi'] in dois_in_file:
                    continue
                item = {
                    'journal_name': journal,
                    'doi': pub['doi'],
                    'title': pub['title'],
                    'link': pub['link'],
                    'published_date': pub['published_date'],
                    'crossref': Journal.get_publication_crossref(pub['doi']),
                }
                publications.append(item)
            if verbose:
                print(f"  Found {len(pubs_info)} publications, {len(publications)} new.")
        Journal.add_publication_to_json(publications, "data/publications.json")
        return publications


    

########################################################################################
# Testing

# Journal.fetch_publications_from_feeds()
# print(Journal.get_dois_from_rss_rsc("http://feeds.rsc.org/rss/sc"))
# dois = Journal.get_dois_from_rss_rsc("http://feeds.rsc.org/rss/CS")
# for doi in dois:
#     pub = Journal.get_publication_crossref(doi)
#     if pub:
#         print(pub['title'], pub['published_date'])

# print(Journal.get_dois_from_rss("http://pubs.acs.org/action/showFeed?jc=jacsat&type=etoc&feed=rss"))
# print(Journal.get_pubmed_ids("J Am Chem Soc[journal]", retmax=50))        
# print(Journal.get_publication_crossref("10.1021/jacs.3c12345"))
# print(Journal.get_publications_pubmed(["38157456", "38157123"]))    


# print(Journal.get_new_pubs_from_rss(Journal.JOURNAL_FEEDS["Advanced Materials"]))
# print(Journal.get_new_pubs_from_rss(Journal.JOURNAL_FEEDS["Angewandte Chemie Int. Ed."]))

# for journal, feed_url in Journal.JOURNAL_FEEDS.items():
#     print(f"Fetching DOIs from {journal}...")
#     pubs_info = Journal.get_new_pubs_from_rss(feed_url)
#     for pub in pubs_info[:10]:
#         print(f"  {pub['doi']} - {pub['title']} ({pub['published_date']}) - {pub['link']}")