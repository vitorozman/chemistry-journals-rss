import streamlit as st
import json
from datetime import datetime, timedelta
from collections import defaultdict


def refresh_feeds():
    """Refresh the RSS feeds and update the publications data."""
    from journals import Journal
    Journal.fetch_publications_from_feeds(verbose=False)
    st.success("Feeds refreshed successfully!")

def load_publications():
    """Load publications from JSON file."""
    try:
        with open("data/publications.json", "r") as f:
            return json.load(f)
    except FileNotFoundError:
        st.error("Publications data file not found!")
        return []

def filter_recent_publications(publications, days=5):
    """Filter publications from the last N days."""
    cutoff_date = datetime.now() - timedelta(days=days)
    recent_pubs = []
    
    for pub in publications:
        if pub.get('published_date'):
            try:
                pub_date = datetime.strptime(pub['published_date'], '%Y-%m-%d')
                if pub_date >= cutoff_date:
                    recent_pubs.append(pub)
            except (ValueError, TypeError):
                # Skip publications with invalid dates
                continue
    
    return recent_pubs

def format_date(date_str):
    """Format date from YYYY-MM-DD to 'D. MMM YYYY'."""
    try:
        date_obj = datetime.strptime(date_str, '%Y-%m-%d')
        return date_obj.strftime('%-d. %b %Y')
    except (ValueError, TypeError):
        return date_str

def group_and_sort_publications(publications):
    """Group publications by journal and sort them."""
    journal_groups = defaultdict(list)
    
    # Group by journal
    for pub in publications:
        journal_name = pub.get('journal_name', 'Unknown Journal')
        journal_groups[journal_name].append(pub)
    
    # Sort journals alphabetically
    sorted_journals = dict(sorted(journal_groups.items()))
    
    # Sort publications within each journal by date (newest first) then title
    for journal in sorted_journals:
        sorted_journals[journal].sort(
            key=lambda x: (
                datetime.strptime(x.get('published_date', '1900-01-01'), '%Y-%m-%d'),
                x.get('title', '').lower()
            ),
            reverse=True
        )
    
    return sorted_journals

def display_publication(pub):
    """Display a single publication in one compact line."""
    # Format the date
    formatted_date = format_date(pub.get('published_date', ''))
    
    # Get publication details
    title = pub.get('title', 'No Title')
    doi = pub.get('doi', '')
    crossref_access = pub.get('crossref', {}).get('access', False)
    
    # Create components for the single line
    date_part = f"**{formatted_date}**"
    title_part = title
    
    
    # DOI link
    doi_part = f"DOI: [{doi}](https://doi.org/{doi})" if doi else ""
    
    # Crossref indicator
    crossref_part = "✓ Crossref" if crossref_access else ""
    
    # Combine everything in one line
    parts = [date_part, title_part, doi_part, crossref_part]
    # Filter out empty parts
    parts = [part for part in parts if part]
    
    # Join with separators
    line = f"{parts[0]} | {parts[1]} | {' | '.join(parts[2:])}"
    st.markdown(line)

# Main app
st.set_page_config(page_title="What's New in Chemistry", page_icon="🧪", layout="wide")

# Add custom CSS for margins and compact layout
st.markdown("""
<style>
.main .block-container {
    padding-left: 5rem;
    padding-right: 5rem;
    max-width: none;
}
.stMarkdown {
    margin-bottom: 0.2rem;
}
h1, h2, h3 {
    margin-bottom: 0.5rem !important;
    margin-top: 1rem !important;
}
</style>
""", unsafe_allow_html=True)

st.title("What's New - Chemistry Journal Publications")
st.markdown("*Publications from the last 2 weeks*")

if st.button("Refresh", key="refresh_feeds"):
    with st.spinner("Refreshing feeds... Please wait."):
        refresh_feeds()
    st.rerun()

# Load and filter publications
publications = load_publications()

if not publications:
    st.warning("No publications data available.")
    st.stop()

recent_publications = filter_recent_publications(publications, days=14)

if not recent_publications:
    st.info("No publications found in the last 2 weeks.")
    st.stop()

# Group and sort publications
journal_groups = group_and_sort_publications(recent_publications)

# Display statistics
total_pubs = len(recent_publications)
total_journals = len(journal_groups)
st.markdown(f"**Found {total_pubs} publications from {total_journals} journals**")
st.markdown("---")

# Display publications by journal
for journal_name in sorted(journal_groups.keys()):
    publications_list = journal_groups[journal_name]
    
    # Journal header
    st.subheader(f"📖 {journal_name}")
    st.markdown(f"*{len(publications_list)} publication{'s' if len(publications_list) != 1 else ''}*")
    
    # Display publications
    for pub in publications_list:
        display_publication(pub)
    
    st.markdown("")  # Add minimal space between journals