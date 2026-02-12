import feedparser
import requests
import time
import json
from datetime import datetime
from Bio import Entrez


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
        "Nature Energy": "https://www.nature.com/nenergy.rss",
        "ACS Nano": "http://pubs.acs.org/action/showFeed?jc=ancac3&type=etoc&feed=rss",
        "ACS Central Science": "https://pubs.acs.org/action/showFeed?jc=acscii&type=etoc&feed=rss",
        "ACS Catalysis": "https://pubs.acs.org/action/showFeed?jc=accacs&type=etoc&feed=rss",
        "Organic Letters": "http://pubs.acs.org/action/showFeed?jc=orlef7&type=etoc&feed=rss",
        "Inorganic Chemistry": "http://pubs.acs.org/action/showFeed?jc=inocaj&type=etoc&feed=rss",
        "Analytical Chemistry": "http://pubs.acs.org/action/showFeed?jc=ancham&type=etoc&feed=rss",
        "J. Phys. Chem. Letters": "http://pubs.acs.org/action/showFeed?jc=jpclcd&type=etoc&feed=rss",
        "Chemistry - A European Journal": "https://onlinelibrary.wiley.com/feed/15213765/most-recent",
        "RCS Analyst": "http://feeds.rsc.org/rss/an",
        "RCS Analytical Methods": "http://feeds.rsc.org/rss/ay",
        "RCS Catalysis Science & Technology": "http://feeds.rsc.org/rss/cy",
        "RCS Chemical Science": "http://feeds.rsc.org/rss/sc",
        "RCS Chemical Society Reviews": "http://feeds.rsc.org/rss/CS",
        "RCS CrystEngComm": "http://feeds.rsc.org/rss/CE",
        "RCS Digital Discovery": "http://feeds.rsc.org/rss/dd",
        "RCS Sustainable Food Technology": "http://feeds.rsc.org/rss/fb",
        "RCS Sustainable Energy & Fuels": "http://feeds.rsc.org/rss/se",
        "RCS Soft Matter": "http://feeds.rsc.org/rss/SM",
        "RCS Sensors & Diagnostics": "http://feeds.rsc.org/rss/sd",
        "RSC Sustainability": "http://feeds.rsc.org/rss/su",
        "RSC Pharmaceutics": "http://feeds.rsc.org/rss/pm",
        "RSC Medicinal Chemistry": "http://feeds.rsc.org/rss/md",
        "RSC Chemical Biology": "http://feeds.rsc.org/rss/cb",
        "RSC Applied Polymers": "http://feeds.rsc.org/rss/lp",
        "RSC Applied Interfaces": "http://feeds.rsc.org/rss/lf",
        "RSC Advances": "http://feeds.rsc.org/rss/ra",
        "RSC Reaction Chemistry & Engineering": "http://feeds.rsc.org/rss/re",
        "RSC Polymer Chemistry": "http://feeds.rsc.org/rss/py",
        "RSC Physical Chemistry Chemical Physics": "http://feeds.rsc.org/rss/CP",
        "RSC Organic Chemistry Frontiers": "http://feeds.rsc.org/rss/qo",
        "RSC Organic & Biomolecular Chemistry": "http://feeds.rsc.org/rss/OB",
        "RSC New Journal of Chemistry": "http://feeds.rsc.org/rss/NJ",
        "RSC Natural Product Reports": "http://feeds.rsc.org/rss/NP",
        "RSC Nanoscale Advances": "http://feeds.rsc.org/rss/na",
        "RSC Nanoscale": "http://feeds.rsc.org/rss/NR",
        "RSC Molecular Systems Design & Engineering": "http://feeds.rsc.org/rss/me",
        "RSC Materials Horizons": "http://feeds.rsc.org/rss/mh",
        "RSC Materials Chemistry Frontiers": "http://feeds.rsc.org/rss/qm",
        "RSC Materials Advances": "http://feeds.rsc.org/rss/ma",
        "RSC Journal of Materials Chemistry C": "http://feeds.rsc.org/rss/tc",
        "RSC Journal of Materials Chemistry B": "http://feeds.rsc.org/rss/tb",
        "RSC Journal of Materials Chemistry A": "http://feeds.rsc.org/rss/ta",
        "RSC Journal of Analytical Atomic Spectrometry": "http://feeds.rsc.org/rss/ja",
        "RSC Inorganic Chemistry Frontiers": "http://feeds.rsc.org/rss/qi",
        "RSC Industrial Chemistry & Materials": "http://feeds.rsc.org/rss/im",
        "RSC Environmental Science: Water Research & Technology": "http://feeds.rsc.org/rss/ew",
        "RSC Environmental Science: Processes & Impacts": "http://feeds.rsc.org/rss/em",
        "RSC Environmental Science: Nano": "http://feeds.rsc.org/rss/en",
        "RSC Environmental Science: Atmospheres": "http://feeds.rsc.org/rss/ea",
        "RSC Environmental Science: Advances": "http://feeds.rsc.org/rss/va",
        "RSC Energy Advances": "http://feeds.rsc.org/rss/ya",
        "RSC Energy & Environmental Science": "http://feeds.rsc.org/rss/EE",
        "RSC EES Solar": "http://feeds.rsc.org/rss/el",
        "RSC EES Catalysis": "http://feeds.rsc.org/rss/ey",
        "RSC EES Batteries": "http://feeds.rsc.org/rss/eb"
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
        pubs_info = []
        for entry in feed.get('entries', []):
            doi = parse_doi_from_summary(entry.get("summary", ""))
            if doi:
                pubs_info.append({
                    'doi': doi,
                    'title': entry.get("title"),
                    'link': entry.get("link"),
                    'published_date': Journal.parse_date_from_feed(entry)
                })
        return pubs_info

    @staticmethod
    def get_new_publications(journal, feed_url):
        if "rsc.org" in feed_url:
            return Journal.get_dois_from_rss_rsc(feed_url)
        else:
            return Journal.get_new_pubs_from_rss(feed_url)


    @staticmethod
    def get_abstract(doi, data):
        """Extract abstract from Crossref metadata if available, otherwise try PubMed."""
        crossref_abstract = data.get('abstract')
        if crossref_abstract:
            return crossref_abstract.strip().replace('\n', ' ')
        pub = Journal.get_publication_pubmed_from_doi(doi)

        if pub.get('access'):
            return pub.get('abstract').strip().replace('\n', ' ')
        return None

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
            'published_date': data.get('published', {}).get('date-parts', [[]])[0],
            'authors': [{
                'first_name': a.get("given", ""),
                'last_name': a.get("family", ""),
                'institution': [aff.get('name') for aff in a.get("affiliation", []) if 'name' in aff],
                'orcid': a.get("ORCID", "").replace("https://orcid.org/", "") if a.get("ORCID") else None
            } for a in data.get("author", [])],
            'abstract': Journal.get_abstract(doi, data)
        }
    

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
    def get_publication_pubmed_from_doi(doi):
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
    def fetch_publications_from_feeds(verbose=True):
        dois_in_file = set(p.get('doi') for p in Journal.load_from_json("data/publications.json") or [])
        for journal, feed_url in Journal.JOURNAL_FEEDS.items():
            publications = []
            if verbose:
                print(f"Fetching publications from {journal}...")
            pubs_info = Journal.get_new_publications(journal, feed_url)
            # pubs_info = Journal.get_new_pubs_from_rss(feed_url)
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




# print(Journal.get_publication_crossref('10.1073/pnas.2413938121'))
# print(Journal.get_publication_crossref('10.1021/acs.accounts.5c00783'))


# Journal.fetch_publications_from_feeds()

# for journal, feed_url in Journal.JOURNAL_FEEDS_EXTENDED.items():
#     print(f"Fetching DOIs from {journal}...")
#     pubs_info = Journal.get_new_publications(journal, feed_url)
#     print(f"  Found {len(pubs_info)} publications.")
#     if len(pubs_info) == 0:
#         print(f"    Check {journal}: link {feed_url}")
#     else: 
#         for pub in pubs_info[:1]:
#             print(f"  {pub['doi']} - {pub['title']} ({pub['published_date']}) - {pub['link']}")

# print(Journal.get_dois_from_rss_rsc("http://feeds.rsc.org/rss/sc"))
# for journal, feed_url in Journal.JOURNAL_FEEDS.items():
#     print(f"Fetching DOIs from {journal}...")
#     if "rsc.org" in feed_url:
#         pubs = Journal.get_dois_from_rss_rsc(feed_url)
#         print(f"  Found {len(pubs)} publications.")
#         if len(pubs) == 0:
#             print(f"    Check {journal}: link {feed_url}")
#     else:
#         pubs = Journal.get_new_pubs_from_rss(feed_url) 
#         print(f"  Found {len(pubs)} publications.")
#         if len(pubs) == 0:
#             print(f"    Check {journal}: link {feed_url}")


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