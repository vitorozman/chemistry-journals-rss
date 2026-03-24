import re
from datetime import datetime

import feedparser
import requests


class BaseRSSParser:
    """Generic RSS/Atom parser for journal feeds."""

    DOI_PATTERN = re.compile(r"10\.\d{4,9}/[-._;()/:A-Z0-9]+", re.IGNORECASE)
    DATE_FORMATS = (
        "%a, %d %b %Y %H:%M:%S %z",
        "%Y-%m-%dT%H:%M:%SZ",
        "%Y-%m-%d",
    )
    def parse(self, feed_url, verify_ssl=True):
        feed = self._fetch_feed(feed_url, verify_ssl=verify_ssl)
        if not feed or "entries" not in feed:
            return []

        pubs_info = []
        for entry in feed.get("entries", []):
            pubs_info.append(
                {
                    "doi": self.extract_doi(entry),
                    "title": self._clean_text(entry.get("title")),
                    "link": self._clean_text(entry.get("link")),
                    "published_date": self.extract_date(entry),
                }
            )
        return pubs_info

    def extract_doi(self, entry):
        for key in ("prism_doi", "dc_identifier", "doi"):
            doi = self.normalize_doi(entry.get(key))
            if doi:
                return doi

        link = entry.get("link", "")
        # print(self.normalize_doi(link))
        return self.normalize_doi(link)

    def extract_date(self, entry):
        for key in ("prism_coverdate", "dc_date", "updated", "published"):
            parsed = self._parse_date_text(entry.get(key))
            if parsed:
                return parsed

        for key in ("published_parsed", "updated_parsed"):
            date_struct = entry.get(key)
            if date_struct:
                return f"{date_struct.tm_year:04d}-{date_struct.tm_mon:02d}-{date_struct.tm_mday:02d}"

        return None

    def _fetch_feed(self, feed_url, verify_ssl=True):
        feed = None
        try:
            response = requests.get(
                feed_url,
                timeout=10,
                headers={"User-Agent": "Mozilla/5.0"},
                verify=verify_ssl,
            )
            response.raise_for_status()
            feed = feedparser.parse(response.content)
        except requests.exceptions.SSLError:
            if verify_ssl:
                try:
                    response = requests.get(
                        feed_url,
                        timeout=10,
                        headers={"User-Agent": "Mozilla/5.0"},
                        verify=False,
                    )
                    response.raise_for_status()
                    feed = feedparser.parse(response.content)
                except Exception:
                    feed = None
        except Exception:
            feed = None

        if feed is None:
            try:
                feed = feedparser.parse(feed_url)
            except Exception:
                feed = None

        return feed

    def _parse_date_text(self, date_text):
        date_text = self._clean_text(date_text)
        if not date_text:
            return None

        if len(date_text) >= 10:
            first_10 = date_text[:10]
            if re.fullmatch(r"\d{4}-\d{2}-\d{2}", first_10):
                return first_10

        for fmt in self.DATE_FORMATS:
            try:
                return datetime.strptime(date_text, fmt).strftime("%Y-%m-%d")
            except ValueError:
                continue

        return None

    @classmethod
    def normalize_doi(cls, raw_value):
        value = cls._clean_text(raw_value)
        if not value:
            return None

        value = re.sub(r"^doi\s*:\s*", "", value, flags=re.IGNORECASE)
        value = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", value, flags=re.IGNORECASE)

        doi_match = cls.DOI_PATTERN.search(value)
        if not doi_match:
            return None

        return doi_match.group(0).rstrip(".,;)")

    @staticmethod
    def _clean_text(value):
        if value is None:
            return None
        if isinstance(value, str):
            return value.strip()
        return str(value).strip()


class RSCRSSParser(BaseRSSParser):
    """RSC feeds usually place DOI inside summary HTML."""

    RSC_SUMMARY_PATTERN = re.compile(
        r"(?:DOI|<b>DOI</b>)\s*:\s*(10\.\d{4,9}/[^\s,<]+)", re.IGNORECASE
    )

    def extract_doi(self, entry):
        summary = self._clean_text(entry.get("summary")) or ""
        match = self.RSC_SUMMARY_PATTERN.search(summary)
        if match:
            return match.group(1)
        return super().extract_doi(entry)


class AnnualReviewRSSParser(BaseRSSParser):
    """Annual Reviews feeds often require DOI extraction from multiple identifier fields."""

    def extract_doi(self, entry):
        for key in ("dc_identifier", "id", "link", "summary"):
            doi = self.normalize_doi(entry.get(key))
            if doi:
                return doi
        return super().extract_doi(entry)


class ELifeRSSParser(BaseRSSParser):
    """eLife feeds typically include DOI only in the id field."""

    def extract_doi(self, entry):
        doi = self.normalize_doi(entry.get("id"))
        if doi:
            return doi
        return super().extract_doi(entry)


class RSSParserFactory:
    """Returns the right parser implementation for each journal/feed."""

    _default_parser = BaseRSSParser()
    _rsc_parser = RSCRSSParser()
    _annual_review_parser = AnnualReviewRSSParser()
    _elife_parser = ELifeRSSParser()

    _journal_parsers = {
        "Annual Review of Analytical Chemistry": _annual_review_parser,
        "Annual Review of Biochemistry": _annual_review_parser,
        "eLife": _elife_parser,
    }

    @classmethod
    def register_journal_parser(cls, journal_name, parser):
        cls._journal_parsers[journal_name] = parser

    @classmethod
    def get_parser(cls, journal_name, feed_url):
        if journal_name in cls._journal_parsers:
            return cls._journal_parsers[journal_name]
        if "rsc.org" in feed_url:
            return cls._rsc_parser
        return cls._default_parser


if __name__ == "__main__":
    # Example usage: parser is selected by journal name/feed URL.
    example_feeds = {
        "Chemical Engineering Journal": "https://rss.sciencedirect.com/publication/science/13858947",
        # "Annual Review of Analytical Chemistry": "https://www.annualreviews.org/action/showFeed?ui=45mu4&mi=3fndc3&ai=6690&jc=anchem&type=etoc&feed=atom",
        # "Nature": "https://www.nature.com/nature.rss",
    }

    for journal_name, feed_url in example_feeds.items():
        parser = RSSParserFactory.get_parser(journal_name, feed_url)
        publications = parser.parse(feed_url)

        print(
            f"{journal_name} -> {parser.__class__.__name__} -> {len(publications)} entries"
        )
        if publications:
            print(publications[0])
