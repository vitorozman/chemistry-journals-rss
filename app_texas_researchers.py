import streamlit as st
import json
from collections import defaultdict

PAGE_SIZE = 100  # researchers per page


# ── Data loading ──────────────────────────────────────────────────────────────

@st.cache_data
def load_candidates() -> list[dict]:
    """
    Load Texas researcher candidates and pre-compute lowercase search fields
    so filtering never rebuilds strings on every query.
    """
    try:
        with open("data/texas_candidates.json", "r") as f:
            raw = json.load(f)
    except FileNotFoundError:
        st.error("texas_candidates.json not found in data/")
        return []

    indexed = []
    for r in raw:
        indexed.append({
            **r,
            "_first": r.get("first_name", "").lower(),
            "_last": r.get("last_name", "").lower(),
            "_insts": r.get("primary_institution", "").lower(),
        })
    return indexed


# ── Filtering ─────────────────────────────────────────────────────────────────

def filter_candidates(indexed: list[dict], query: str) -> list[dict]:
    """
    Filter using pre-computed search fields. Each token must appear as a
    case-insensitive substring in first name, last name, or any institution.
    """
    if not query:
        return indexed

    tokens = query.lower().split()
    return [
        r for r in indexed
        if all(
            t in r["_first"] or t in r["_last"] or t in r["_insts"]
            for t in tokens
        )
    ]


# ── Grouping / sorting ────────────────────────────────────────────────────────

def group_by_institution(candidates: list[dict]) -> dict[str, list[dict]]:
    """Group researchers by primary_institution, alphabetically sorted."""
    groups: dict[str, list] = defaultdict(list)
    for r in candidates:
        inst = r.get("primary_institution") or "Unknown Institution"
        groups[inst].append(r)
    for inst in groups:
        groups[inst].sort(key=lambda r: (r.get("last_name", "").lower(),
                                         r.get("first_name", "").lower()))
    return dict(sorted(groups.items()))


# ── HTML builders (no Streamlit widgets per researcher) ───────────────────────

def _pub_html(pub: dict) -> str:
    title = pub.get("title", "Untitled").replace("\u00a0", " ")
    doi = pub.get("doi", "")
    year = pub.get("publication_date", "").split("-")[0]
    title_html = (
        f'<a href="https://doi.org/{doi}" target="_blank">{title}</a>' if doi else title
    )
    year_html = f' <span class="pub-date">{year}</span>' if year else ""
    return f"<li>{title_html}{year_html}</li>"


def _researcher_html(r: dict) -> str:
    first = r.get("first_name", "")
    last = r.get("last_name", "")
    full_name = f"{first} {last}".strip()
    pubs = r.get("publications", [])
    pub_count = len(pubs)
    pub_items = "".join(_pub_html(p) for p in pubs)
    pub_list = f"<ol class='pub-list'>{pub_items}</ol>" if pubs else "<p class='no-pubs'>No publications.</p>"
    return (
        f'<details>'
        f'<summary><span class="researcher-name">{full_name}</span>'
        f'<span class="pub-badge">Pubs {pub_count}</span></summary>'
        f'<div class="pub-body">{pub_list}</div>'
        f'</details>'
    )


def render_institution_block(institution: str, researchers: list[dict]) -> None:
    """Render a whole institution as ONE st.markdown call."""
    count = len(researchers)
    header = (
        f'<div class="inst-header">'
        f'<span class="inst-name">{institution}</span>'
        f'<span class="inst-count">{count} researcher{"s" if count != 1 else ""}</span>'
        f'</div>'
    )
    researchers_html = "".join(_researcher_html(r) for r in researchers)
    st.markdown(
        f'{header}<div class="inst-body">{researchers_html}</div>'
        f'<div class="inst-divider"></div>',
        unsafe_allow_html=True,
    )


# ── CSS ───────────────────────────────────────────────────────────────────────

CUSTOM_CSS = """
<style>
/* ── page layout ── */
.main .block-container {
    padding: 2rem 4rem 3rem 4rem;
    max-width: 1100px;
}

/* ── search bar ── */
div[data-testid="stTextInput"] input {
    font-size: 1rem;
    border-radius: 8px;
}

/* ── institution header ── */
.inst-header {
    display: flex;
    align-items: baseline;
    gap: 0.6rem;
    margin-top: 1.8rem;
    margin-bottom: 0.3rem;
}
.inst-name {
    font-size: 1.05rem;
    font-weight: 700;
    color: #111827;
    letter-spacing: 0.01em;
}
.inst-count {
    font-size: 0.75rem;
    font-weight: 400;
    color: #6b7280;
}
.inst-body {
    margin-left: 0.5rem;
}
.inst-divider {
    border-top: 1px solid #e5e7eb;
    margin-top: 0.8rem;
}

/* ── researcher row (native <details>) ── */
details {
    border-bottom: 1px solid #f3f4f6;
    padding: 0.15rem 0;
}
details summary {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    cursor: pointer;
    list-style: none;
    padding: 0.2rem 0;
    user-select: none;
}
details summary::-webkit-details-marker { display: none; }
details summary::before {
    content: "▶";
    font-size: 0.6rem;
    color: #9ca3af;
    flex-shrink: 0;
    transition: transform 0.15s;
}
details[open] summary::before {
    transform: rotate(90deg);
}
.researcher-name {
    font-size: 0.92rem;
    font-weight: 600;
    color: #1f2937;
    flex: 1;
}
details summary:hover .researcher-name {
    color: #2563eb;
}
.pub-badge {
    font-size: 0.7rem;
    font-weight: 500;
    color: #6b7280;
    background: #f3f4f6;
    border: 1px solid #e5e7eb;
    border-radius: 20px;
    padding: 1px 7px;
    flex-shrink: 0;
}
.pub-body {
    padding: 0.1rem 0 0.3rem 1.2rem;
}

/* ── publication list ── */
ol.pub-list {
    margin: 0;
    padding-left: 1.3rem;
}
ol.pub-list li {
    font-size: 0.82rem;
    color: #374151;
    line-height: 1.45;
    margin-bottom: 0.12rem;
}
ol.pub-list li a {
    color: #2563eb;
    text-decoration: none;
}
ol.pub-list li a:hover { text-decoration: underline; }
.pub-date {
    color: #9ca3af;
    font-size: 0.75rem;
    margin-left: 4px;
}
.no-pubs {
    font-size: 0.8rem;
    color: #9ca3af;
    font-style: italic;
    margin: 0;
}

/* ── summary bar ── */
.summary-bar {
    font-size: 0.85rem;
    color: #6b7280;
    margin-bottom: 0.5rem;
}

/* ── hide Streamlit branding ── */
#MainMenu, footer { visibility: hidden; }
</style>
"""


# ── App entry point ───────────────────────────────────────────────────────────

def main() -> None:
    st.set_page_config(
        page_title="Texas Researchers Candidates",
        page_icon="🔬",
        layout="centered",
    )
    st.markdown(CUSTOM_CSS, unsafe_allow_html=True)

    # ── Title
    st.title("Texas Researchers Candidates")

    # ── Global search bar
    query = st.text_input(
        label="Search",
        placeholder="🔍  Search by name or institution…",
        label_visibility="collapsed",
        key="search_query",
    )

    # ── Load data (cached, includes pre-computed search index)
    indexed = load_candidates()
    if not indexed:
        st.stop()

    # ── Filter (fast: pre-computed fields) + group
    filtered = filter_candidates(indexed, query.strip())
    grouped = group_by_institution(filtered)

    # ── Summary
    total_researchers = len(filtered)
    total_insts = len(grouped)
    st.markdown(
        f'<div class="summary-bar">'
        f'<b>{total_researchers}</b> researcher{"s" if total_researchers != 1 else ""} · '
        f'<b>{total_insts}</b> institution{"s" if total_insts != 1 else ""}'
        f'<br><span style="font-size:0.8rem;color:#9ca3af;">Texas researchers found in author list of current researchers in db which could be possible candidates for the investigation</span>'
        f'</div>',
        unsafe_allow_html=True,
    )

    if not filtered:
        st.info("No researchers match your search.")
        st.stop()

    st.markdown("---")

    # ── Pagination over institutions
    inst_names = list(grouped.keys())
    total_pages = max(1, -(-len(inst_names) // PAGE_SIZE))  # ceiling division

    if "page" not in st.session_state:
        st.session_state.page = 0
    # Reset to page 0 whenever the query changes
    if st.session_state.get("_last_query") != query:
        st.session_state.page = 0
        st.session_state["_last_query"] = query

    page = st.session_state.page
    page_insts = inst_names[page * PAGE_SIZE : (page + 1) * PAGE_SIZE]

    # ── Render current page — one st.markdown call per institution
    for inst in page_insts:
        render_institution_block(inst, grouped[inst])

    # ── Pagination controls (only shown when needed)
    if total_pages > 1:
        st.markdown("---")
        col_prev, col_info, col_next = st.columns([1, 4, 1])
        with col_prev:
            if st.button("← Prev", disabled=page == 0):
                st.session_state.page -= 1
                st.rerun()
        with col_info:
            st.markdown(
                f'<div style="text-align:center;font-size:0.85rem;color:#6b7280;padding-top:0.4rem">'
                f'Page {page + 1} of {total_pages} &nbsp;·&nbsp; '
                f'institutions {page * PAGE_SIZE + 1}–{min((page + 1) * PAGE_SIZE, total_insts)} of {total_insts}'
                f'</div>',
                unsafe_allow_html=True,
            )
        with col_next:
            if st.button("Next →", disabled=page >= total_pages - 1):
                st.session_state.page += 1
                st.rerun()


if __name__ == "__main__":
    main()
