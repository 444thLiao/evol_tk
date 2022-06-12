"""
Implement the function for mergeing multiple pdf file into a single page

Currently, it only handles 4/6 figs at total. And the input figs should be the plotly.graph_objs instead of pdf files, directly.

More pdf files into a single page should be carefully handle the layout.


It mainly used in `taxon specific.ipynb of machine learning project` currently.

"""

from PyPDF2 import PdfFileReader, PdfFileWriter, PdfFileMerger
from PyPDF2.pdf import PageObject
from decimal import Decimal
import io,os
from pdf2image import convert_from_path

def out_four(figs, ofile=None,figsize=200):
    figs_pdf = []
    for fig in figs:
        r = PdfFileReader(io.BytesIO(fig.to_image('pdf')))
        figs_pdf.append(r.getPage(0))
    if len(figs_pdf) == 4:
        w = figs_pdf[0].mediaBox.getWidth()
        h = figs_pdf[0].mediaBox.getHeight()
        aw = 2*w - Decimal(0.2)*w
        ah = 2*h - Decimal(0.2)*h
        translated_page = PageObject.createBlankPage(None, aw, ah)
        # mergeScaledTranslatedPage(page2, rotation, tx, ty, expand=True)
        translated_page.mergeScaledTranslatedPage(figs_pdf[0], 1, 0, Decimal(0.8)*h)  
        translated_page.mergeScaledTranslatedPage(figs_pdf[1], 1, Decimal(0.8)*w, Decimal(0.8)*h)  
        translated_page.mergeScaledTranslatedPage(figs_pdf[2], 1, 0, 0)  
        translated_page.mergeScaledTranslatedPage(figs_pdf[3], 1, Decimal(0.8)*w, 0)  
    elif len(figs_pdf) == 6:
        w = figs_pdf[0].mediaBox.getWidth()
        h = figs_pdf[0].mediaBox.getHeight()
        aw = 2*w - Decimal(0.2)*w
        ah = 3*h - Decimal(0.4)*h
        translated_page = PageObject.createBlankPage(None, aw, ah)
        # mergeScaledTranslatedPage(page2, rotation, tx, ty, expand=True)
        translated_page.mergeScaledTranslatedPage(figs_pdf[0], 1, 0, Decimal(1.6)*h)  
        translated_page.mergeScaledTranslatedPage(figs_pdf[1], 1, Decimal(0.8)*w, Decimal(1.6)*h)  
        translated_page.mergeScaledTranslatedPage(figs_pdf[2], 1, 0, Decimal(0.8)*h)  
        translated_page.mergeScaledTranslatedPage(figs_pdf[3], 1, Decimal(0.8)*w, Decimal(0.8)*h)  
        translated_page.mergeScaledTranslatedPage(figs_pdf[4], 1, 0, 0)  
        translated_page.mergeScaledTranslatedPage(figs_pdf[5], 1, Decimal(0.8)*w, 0)     
    pdf_write = PdfFileWriter()
    pdf_write.addPage(translated_page)
    if ofile is None:
        ofile = './tmp.pdf'
        with open(ofile, 'wb') as fh:
            pdf_write.write(fh)
        pages = convert_from_path(ofile, figsize)
        os.system(f"rm {ofile}")
        pages[0].save('./tmp.png','PNG')
        
        bytes_fig = open('./tmp.png','rb').read()
        os.system(f"rm ./tmp.png")
        return bytes_fig
    else:
        with open(ofile, 'wb') as fh:
            pdf_write.write(fh)