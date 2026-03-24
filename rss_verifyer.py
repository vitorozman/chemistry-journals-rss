
import feedparser
import requests







def verify_rss_feed(url, verify_ssl=True):

    feed = None
    try:
        response = requests.get(url, timeout=10, headers={'User-Agent': 'Mozilla/5.0'}, verify=verify_ssl)
        response.raise_for_status()
        feed = feedparser.parse(response.content)
    except requests.exceptions.SSLError:
        if verify_ssl:
            try:
                response = requests.get(url, timeout=10, headers={'User-Agent': 'Mozilla/5.0'}, verify=False)
                response.raise_for_status()
                feed = feedparser.parse(response.content)
            except:
                feed = None
    except:
        feed = None
    try:
        feed = feedparser.parse(url)
    except:
        feed = None
    
    if not feed or 'entries' not in feed:
        return []
    pubs = []
    for entry in feed.get('entries', []):
        title = entry.get('title', 'No Title')
        pubs.append(title)
    return pubs



########### Example usage:

JOURNALS_FEED_NEW = {
    
    "Nature": "https://www.nature.com/nature.rss",  
    
}

if __name__ == "__main__":
    for journal, url in JOURNALS_FEED_NEW.items():
        print(f"Verifying RSS feed for {journal}...")
        publications = verify_rss_feed(url)
        if publications:
            print(f"Publications found in RSS feed: {len(publications)}")
            for pub in publications[:5]:  # Print first 5 publications
                print(f"- {pub}")
        else:
            print("No publications found or failed to parse RSS feed.")
        print("\n")
        input("Press Enter to continue to the next journal...")


