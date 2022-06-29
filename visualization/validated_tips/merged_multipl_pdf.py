"""
Implement the function for mergeing multiple pdf file into a single page

Currently, it only handles 4/6 figs at total. And the input figs should be the plotly.graph_objs instead of pdf files, directly.

More pdf files into a single page should be carefully handle the layout.


It mainly used in `taxon specific.ipynb of machine learning project` currently.

The second section of scripts are used in anammox project

"""

from PyPDF2 import PdfFileReader, PdfFileWriter, PdfFileMerger
from PyPDF2.pdf import PageObject
from decimal import Decimal
import io,os
from pdf2image import convert_from_path
from tqdm import tqdm
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics

pdfmetrics.registerFont(TTFont('Arial-ita', '/home-user/thliao/.fonts/Arial Italic.ttf'))

def merged_pdfs(pdf_list,ofile='./tmp.pdf' ,shape=3):
    """
    It could only merged pdf into image. like png,jpeg
    
     shape is either a tuple or a intergrate. If it is a intergrate, it is taken as the number of cols. The number of cols will be calculated according to the number of pdf_list
     if it is a tuple, it shoule be (nrow,ncol)
    """
    if type(shape) == int:
        suffix = 1 if len(pdf_list)%shape else 0
        nrow = len(pdf_list)//shape + suffix
        ncol = shape
    else:
        nrow,ncol = shape
    print(nrow,ncol)
    figs_pdf = []
    names = []
    for pdf in pdf_list:
        r = PdfFileReader(open(pdf,'rb'))
        figs_pdf.append(r.getPage(0))
        names.append(pdf.split('/')[-1].split('_')[0])
    w = figs_pdf[0].mediaBox.getWidth()
    h = figs_pdf[0].mediaBox.getHeight()
    aw = ncol*w; ah = nrow* (h )
    final_page = PageObject.createBlankPage(None, aw, ah)
    pos = (0,0)  #  x,y ; width, height
    for name,pdf_obj in tqdm(zip(names,figs_pdf)):
        #text_in_pdf = make_pdf_text(name,size=(10,10))
        fig_pos = (pos[0]*w, pos[1]*(h))
        final_page.mergeScaledTranslatedPage(pdf_obj, 1, fig_pos[0],fig_pos[1])  
        # mid_x = (w-10)/2
        # final_page.mergeScaledTranslatedPage(text_in_pdf, 1, fig_pos[0]+mid_x,fig_pos[1]+h+50)  
        if pos[0] == ncol-1:
            pos = (0,pos[1]+1)
        else:
            pos = (pos[0]+1,pos[1]) 
    pdf_write = PdfFileWriter()
    pdf_write.addPage(final_page)
    with open(ofile, 'wb') as fh:
        pdf_write.write(fh)
        
def make_pdf_text(text,size=(10,10),fontsize=12):
    packet = io.BytesIO()
    can = canvas.Canvas(packet,pagesize=size)
    can.setFont("Arial-ita",fontsize)
    can.drawCentredString(size[0]/2,size[1]/2,text)
    can.save()
    packet.seek(0)
    return PdfFileReader(packet).getPage(0) 

def merged_pdfs_with_texts(pdf_list,ofile='./tmp.pdf' ,shape=3,out_name=True):
    """
    It could only merged pdf into image. like png,jpeg
    
     shape is either a tuple or a intergrate. If it is a intergrate, it is taken as the number of cols. The number of cols will be calculated according to the number of pdf_list
     if it is a tuple, it shoule be (nrow,ncol)
    """
    if type(shape) == int:
        suffix = 1 if len(pdf_list)%shape else 0
        nrow = len(pdf_list)//shape + suffix
        ncol = shape
    else:
        nrow,ncol = shape
    #print(nrow,ncol)
    figs_pdf = []
    names = []
    for pdf in pdf_list:
        r = PdfFileReader(open(pdf,'rb'))
        figs_pdf.append(r.getPage(0))
        names.append(pdf.split('/')[-1].split('_')[0])
    w = figs_pdf[0].mediaBox.getWidth()
    h = figs_pdf[0].mediaBox.getHeight()
    text_height = 80
    aw = ncol*w; ah = nrow* (h + text_height )
    final_page = PageObject.createBlankPage(None, aw, ah)
    pos = (0,0)  #  x,y ; width, height
    for name,pdf_obj in tqdm(zip(names,figs_pdf)):
        text_in_pdf = make_pdf_text(name,size=(float(w),text_height),fontsize=50)
        fig_pos = (pos[0]*w, pos[1]*(h) + text_height*(pos[1]+1) )
        final_page.mergeScaledTranslatedPage(pdf_obj, 1, fig_pos[0],fig_pos[1])  
        
        final_page.mergeScaledTranslatedPage(text_in_pdf, 1, fig_pos[0],fig_pos[1]-text_height)  
        if pos[0] == ncol-1:
            pos = (0,pos[1]+1)
        else:
            pos = (pos[0]+1,pos[1]) 
    pdf_write = PdfFileWriter()
    pdf_write.addPage(final_page)
    with open(ofile, 'wb') as fh:
        pdf_write.write(fh)
         
from decimal import Decimal
import io,os
from PIL import ImageDraw,ImageFont,Image
from pdf2image import convert_from_path
os.chdir('gene_GandL/genetrees/iqtree_o_pdf/')

def get_figs(pdf_list,outfile='./tmp.png' ,shape=3,each_size=625,title_height=100):
    """
    It could only merged pdf into image. like png,jpeg
    
     shape is either a tuple or a intergrate. If it is a intergrate, it is taken as the number of cols. The number of cols will be calculated according to the number of pdf_list
     if it is a tuple, it shoule be (nrow,ncol)
    """
    if type(shape) == int:
        suffix = 1 if len(pdf_list)%shape else 0
        nrow = len(pdf_list)//shape + suffix
        ncol = shape
    else:
        nrow,ncol = shape
    print(nrow,ncol)
    init_bg = Image.new('RGB',
                        (each_size*ncol, (each_size+title_height)*nrow),  # width,height
                        '#ffffff')
    draw = ImageDraw.Draw(init_bg)
    fillColor = "#000000"
    font = ImageFont.truetype("/home-user/thliao/.fonts/Arial Italic.ttf", 45)
    pos = (0,0)  #  x,y ; width, height
    for f in tqdm(pdf_list):
        images = convert_from_path(f)
        gene = f.split('/')[-1].split('_')[0]
        a = images[0]
        width,height = a.size
        a = a.resize((width//4,height//4))
        fig_pos = (pos[0] * each_size, pos[1]*(each_size+title_height))
        init_bg.paste(a, box= fig_pos )
        text = f"{gene}"
        w,h = font.getsize(text)
        mid_x = (each_size-w)/2
        position = (fig_pos[0] + mid_x, fig_pos[1]+each_size+15)
        draw.text(position,text,fill=fillColor,font=font,anchor ="lb",align ='center')
        if pos[0] == nrow:
            pos = (0,pos[1]+1)
        else:
            pos = (pos[0]+1,pos[1]) 
    init_bg.save(outfile)


