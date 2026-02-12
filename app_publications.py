import streamlit as st
import json
from datetime import datetime, timedelta



def create_pretty_print(publication):
    """
    Create a single-line APA 7 style citation with bold Texas authors,
    clickable DB authors, and a compact summary for Streamlit.
    """
    crossref = publication.get("crossref", {})
    authors = crossref.get("authors", [])
    
    # Format publication date as "3 Jan 2026"
    published_date = crossref.get("published_date", [])
    if published_date and len(published_date) >= 3:
        year, month, day = published_date[0], published_date[1], published_date[2]
        month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                      "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
        formatted_date = f"{day} {month_names[month-1]} {year}"
    elif published_date and len(published_date) >= 1:
        formatted_date = str(published_date[0])  # Just year
    else:
        # Fallback to string date
        fallback_date = publication.get("published_date", "n.d.")
        if fallback_date and fallback_date != "n.d.":
            try:
                from datetime import datetime
                dt = datetime.strptime(fallback_date[:10], '%Y-%m-%d')
                month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
                formatted_date = f"{dt.day} {month_names[dt.month-1]} {dt.year}"
            except:
                formatted_date = "n.d."
        else:
            formatted_date = "n.d."
    
    # Format authors
    formatted_authors = []
    for author in authors:
        name = f"{author['last_name']} {author['first_name'][0]}."
        # Bold if Texas affiliation
        if author.get("texas_affiliation"):
            name = f"**{name}**"
        # Make clickable if DB researcher
        if author.get("db_researcher_id"):
            researcher_id = author["db_researcher_id"]
            name = f"[{name}](https://welch-frontend.fly.dev/researcher/{researcher_id})"
        formatted_authors.append(name)
    
    authors_str = ", ".join(formatted_authors)
    
    # APA 7 style: Author(s). (Year). Title. Journal, doi
    title = crossref.get("title") or publication.get("title")
    # Clean up HTML tags in title for proper display
    if title and ("<sub>" in title or "<sup>" in title or "<i>" in title or "<b>" in title):
        # Keep the HTML tags as they are for proper chemical formula display
        formatted_title = title
    else:
        formatted_title = title
    
    journal = publication.get("journal_name")
    doi = publication.get("doi")
    link = publication.get("link", f"https://doi.org/{doi}")
    
    citation = f"{authors_str} ({formatted_date}). {formatted_title}. *{journal}*. [DOI]({link})"
    
    # Summary in smaller font (display only if present and not empty)
    summary = publication.get("summary", "")
    
    # Streamlit display with compact spacing
    st.markdown(citation, unsafe_allow_html=True)
    if summary and summary.strip():  # Check for both existence and non-empty content
        st.markdown(f"<sub style='margin-top: -0.5rem; display: block;'>{summary}</sub>", unsafe_allow_html=True)

def load_publications():
    """Load publications from JSON file."""
    try:
        with open("data/publications_display.json", "r") as f:
            return json.load(f)
    except FileNotFoundError:
        st.error("Publications with summaries data file not found!")
        return []


def parse_and_sort_publications(publications):
    """Sort publications from newest to oldest by date only."""
    def get_sort_date(pub):
        crossref = pub.get("crossref", {})
        published_date = crossref.get("published_date", [])
        
        # If crossref published_date exists and has at least year
        if published_date and len(published_date) > 0:
            # Pad with zeros if month/day are missing: [2024] -> [2024, 1, 1]
            year = published_date[0] if len(published_date) > 0 else 0
            month = published_date[1] if len(published_date) > 1 else 1
            day = published_date[2] if len(published_date) > 2 else 1
            return (year, month, day)
        
        # Fallback to string date if available
        fallback_date = pub.get('published_date', '')
        if fallback_date:
            try:
                # Parse YYYY-MM-DD format
                parts = fallback_date.split('-')
                return (int(parts[0]), int(parts[1]), int(parts[2]))
            except (ValueError, IndexError):
                pass
        
        # Default to very old date for missing dates
        return (0, 0, 0)
    
    # Sort publications by date: newest first
    sorted_pubs = sorted(publications, key=get_sort_date, reverse=True)
    
    return sorted_pubs

# Main app
st.set_page_config(page_title="What's New in Chemistry", page_icon="🧪", layout="centered")

# Add custom CSS for margins and compact layout
st.markdown("""
<style>
.main .block-container {
    padding-left: 5rem;
    padding-right: 5rem;
    max-width: none;
}
.stMarkdown {
    margin-bottom: 0.1rem;
    margin-top: 0.1rem;
}
h1, h2, h3 {
    margin-bottom: 0.5rem !important;
    margin-top: 1rem !important;
}
.publication-separator {
    margin: 0.5rem 0 !important;
    border-color: #e0e0e0 !important;
}
sub {
    margin-top: -0.3rem !important;
    line-height: 1.2 !important;
}
</style>
""", unsafe_allow_html=True)

st.title("What's New")


# Load all publications
publications = load_publications()

if not publications:
    st.warning("No publications data available.")
    st.stop()

# Sort all publications from newest to oldest
sorted_publications = parse_and_sort_publications(publications)

# Display statistics
total_pubs = len(sorted_publications)
st.markdown(f"**Authors from Texas authored {total_pubs} publications in the top-rated chemistry journals.**")
st.markdown("---")

# Display all publications
for i, pub in enumerate(sorted_publications, 1):
    create_pretty_print(pub)
    if i < len(sorted_publications):  # Don't add separator after last publication
        st.markdown("<hr class='publication-separator'>", unsafe_allow_html=True)