from setuptools import setup, find_packages
setup(
    name = "DocLinks",
    author = "Manuel Holtgrewe",
    license = "THE BEER-WARE LICENSE (Revision 42)",
    url = "http://www.seqan.de/",
    description = "Trac plugin for creating links to the SeqAn Documentation",
    version = "0.1",
    packages = find_packages(exclude=['*.tests*']),
    entry_points = """
        [trac.plugins]
        doc_links = DocLinks.doc_links
    """,
)
