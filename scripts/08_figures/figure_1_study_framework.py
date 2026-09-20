from pathlib import Path
from xml.sax.saxutils import escape

OUT = Path("results/figure_exports")
OUT.mkdir(parents=True, exist_ok=True)
W, H = 1295, 1214

INK = "#10203B"
BLUE = "#2C7FB8"
RED = "#C44E3F"
AMBER = "#C98720"
PURPLE = "#6F4AA8"
TEAL = "#2A9D8F"
GREY = "#63748A"
LIGHT_BLUE = "#EEF6FB"
LIGHT_RED = "#FFF0EE"
LIGHT_PURPLE = "#F6F0FF"
LIGHT_TEAL = "#ECF8F5"
LIGHT_AMBER = "#FFF8E8"
LINE = "#B7C7D6"

svg = []
def add(s): svg.append(s)
def rect(x, y, w, h, fill="white", stroke=LINE, sw=2, rx=12):
    add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')
def line(x1, y1, x2, y2, stroke=LINE, sw=3, dash=None, marker=None):
    d = f' stroke-dasharray="{dash}"' if dash else ''
    m = f' marker-end="url(#{marker})"' if marker else ''
    add(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke}" stroke-width="{sw}"{d}{m}/>')
def text(x, y, s, size=20, fill=INK, weight="400", anchor="start", italic=False):
    style = 'font-style:italic;' if italic else ''
    add(f'<text x="{x}" y="{y}" font-family="Arial,Helvetica,sans-serif" font-size="{size}px" font-weight="{weight}" fill="{fill}" text-anchor="{anchor}" style="{style}">{escape(s)}</text>')
def block(x, y, lines, size=20, fill=INK, weight="400", leading=None, anchor="start"):
    leading = leading or int(size * 1.18)
    add(f'<text x="{x}" y="{y}" font-family="Arial,Helvetica,sans-serif" font-size="{size}px" font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">')
    for i, s in enumerate(lines):
        dy = 0 if i == 0 else leading
        add(f'<tspan x="{x}" dy="{dy}px">{escape(s)}</tspan>')
    add('</text>')
def pill(x, y, w, h, label, fill, size=20):
    rect(x, y, w, h, fill=fill, stroke=fill, sw=1, rx=18)
    text(x+w/2, y+h*0.68, label, size=size, fill="white", weight="700", anchor="middle")
def dot(x, y, fill, r=11, stroke=None):
    st = f' stroke="{stroke}" stroke-width="2"' if stroke else ''
    add(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}"{st}/>')
def square(x, y, fill, r=9):
    add(f'<rect x="{x-r}" y="{y-r}" width="{2*r}" height="{2*r}" rx="2" fill="{fill}"/>')

add(f'''<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="168.65mm" viewBox="0 0 {W} {H}">
<defs><marker id="arrow-blue" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="{BLUE}"/></marker>
<marker id="arrow-red" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="{RED}"/></marker>
<marker id="arrow-amber" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="{AMBER}"/></marker></defs>
<rect width="{W}" height="{H}" fill="white"/>''')

# Header
rect(26, 18, 1243, 118, fill="#F6F9FC", stroke=LINE, sw=2, rx=14)
text(54, 53, "Study question", 20, GREY, "700")
block(54, 83, ["Does an RA-derived transcriptional programme in OA represent a portable disease-like state,",
               "or context-dependent reuse of inflammatory machinery?"], 24, INK, "700", 28)

# Panel a
text(27, 143, "a", 28, INK, "700")
text(64, 143, "External reference and frozen programme", 25, INK, "700")
rect(28, 160, 355, 205, fill=LIGHT_RED, stroke="#E89A91", sw=2)
pill(47, 176, 119, 38, "RA reference", RED, 19)
text(52, 253, "GSE89408 • synovium", 20, RED, "700")
text(52, 290, "152 RA vs 28 healthy", 20, INK, "700")
text(52, 321, "limma–voom differential ranking", 17, GREY)
for i, c in enumerate(["#C44E3F", "#C44E3F", "#E28F86", "#F0B6AF"]): rect(64+i*72, 337, 55, 26, fill=c, stroke=c, sw=0, rx=0)
line(386, 260, 416, 260, stroke=RED, sw=4, marker="arrow-red")

rect(418, 160, 361, 205, fill="#FFF7F5", stroke="#E89A91", sw=2)
text(443, 204, "Pre-validation freeze", 20, RED, "700")
for i, h in enumerate([37, 56, 73, 46, 68, 84, 57, 74, 45]):
    rect(443+i*22, 279-h, 12, h, fill=RED, stroke=RED, sw=0, rx=0)
for i, h in enumerate([40, 56, 70, 46, 64, 78, 53, 69, 43]):
    rect(590+i*18, 279-h, 12, h, fill=BLUE, stroke=BLUE, sw=0, rx=0)
text(496, 310, "150 RA-up", 19, RED, "700", "middle")
text(646, 310, "150 RA-down", 19, BLUE, "700", "middle")
block(600, 336, ["Equal weights • membership fixed", "No OA-driven reselection or reweighting"], 15, GREY, "400", 21, "middle")
line(783, 260, 814, 260, stroke=BLUE, sw=4, marker="arrow-blue")

rect(816, 160, 453, 205, fill="#F7FBFE", stroke=LINE, sw=2)
text(844, 204, "Frozen continuous projection", 20, INK, "700")
text(844, 247, "log₂(TPM+1) → gene-wise z-score", 17, INK)
text(844, 278, "programme score = mean(available genes)", 17, INK)
rect(841, 300, 401, 62, fill="#F4F7FA", stroke=LINE, sw=1.5, rx=10)
text(1042, 331, "RA projection = RA-up − RA-down", 18, RED, "700", "middle")
text(1042, 358, "No threshold • no RA-like / non-RA-like class", 15, GREY, "400", "middle")

# Panel b
text(27, 408, "b", 28, INK, "700")
text(64, 408, "Projection into seven OA cohorts across heterogeneous study contexts", 25, INK, "700")
rect(28, 426, 505, 230, fill=LIGHT_BLUE, stroke="#81B8D5", sw=2)
pill(48, 443, 168, 38, "Principal validation", BLUE, 18)
text(51, 510, "GSE283079 • raw RNA-seq reconstruction", 20, BLUE, "700")
rect(58, 553, 64, 54, fill="white", stroke=BLUE, sw=2, rx=12); rect(58, 545, 64, 16, fill="white", stroke=BLUE, sw=2, rx=8); line(58, 553, 58, 607, BLUE, 2); line(122, 553, 122, 607, BLUE, 2)
line(126, 580, 165, 580, BLUE, 3, marker="arrow-blue")
rect(168, 540, 128, 79, fill="white", stroke="#9CCAE0", sw=2, rx=8); text(232, 572, "QC", 18, BLUE, "700", "middle"); text(232, 598, "paired reads", 14, GREY, "400", "middle"); text(232, 616, "FastQC", 14, GREY, "400", "middle")
line(300, 580, 339, 580, BLUE, 3, marker="arrow-blue")
rect(342, 540, 116, 79, fill="white", stroke="#9CCAE0", sw=2, rx=8); text(400, 572, "Salmon", 18, BLUE, "700", "middle"); text(400, 601, "GENCODE v47", 14, GREY, "400", "middle")
line(462, 580, 478, 580, BLUE, 3, marker="arrow-blue")
add('<circle cx="494" cy="580" r="36" fill="white" stroke="%s" stroke-width="2"/>' % BLUE); text(494, 575, "36 OA", 17, BLUE, "700", "middle"); text(494, 596, "for validation", 13, GREY, "400", "middle")
text(51, 644, "41 SRA runs", 18, BLUE, "700"); text(163, 644, "+ 5 non-OA samples retained for descriptive purposes", 15, GREY)

rect(552, 426, 717, 230, fill="#F9FBFD", stroke=LINE, sw=2)
text(580, 465, "Cross-cohort OA relationship analysis", 20, INK, "700")
cohorts=[("GSE89408", "n=22 • shared study context"), ("GSE55235", "n=10 • GPL96"), ("GSE55457", "n=10 • GPL96"), ("GSE55584", "n=6 • GPL96"), ("GSE82107", "n=10 • GPL570"), ("GSE206848", "n=7 • GPL570")]
xs=[575, 813, 1050]; ys=[482, 557]
for i,(name,desc) in enumerate(cohorts):
    x=xs[i%3]; y=ys[i//3]; rect(x,y,220,58,fill="white",stroke="#A9C4D8",sw=1.5,rx=10); text(x+15,y+25,name,17,BLUE,"700"); text(x+15,y+47,desc,13,GREY)
rect(1020, 620, 222, 35, fill=LIGHT_BLUE, stroke="#A9C4D8", sw=1.5, rx=10); text(1131, 644, "7 OA cohorts • n=101", 16, BLUE, "700", "middle")

# Panel c
text(27, 702, "c", 28, INK, "700")
text(64, 702, "Test the surrounding immune architecture, not the score alone", 25, INK, "700")
rect(28, 720, 505, 226, fill=LIGHT_PURPLE, stroke="#B7A4DB", sw=2)
text(52, 753, "Prespecified comparator programmes", 19, PURPLE, "700")
labels=[("General inflammation",AMBER), ("APC 10-gene",PURPLE), ("MHC-II",PURPLE), ("Myeloid",TEAL), ("Interferon",TEAL), ("T-cell",GREY), ("B-cell",GREY), ("Fibroblast / ECM",GREY), ("Osteoclast",GREY)]
for i,(lab,col) in enumerate(labels):
    x=49+(i%3)*154; y=773+(i//3)*52; rect(x,y,145,42,fill="white",stroke="#B7C7D6",sw=1.5,rx=12); text(x+72,y+27,lab,14,col,"700","middle")
text(52, 928, "Same scoring rule • equal weighting • coverage reported", 14, GREY)

rect(552, 720, 717, 226, fill="#F9FBFD", stroke=LINE, sw=2)
text(580, 754, "Four complementary tests of portability", 20, INK, "700")
rect(1004, 735, 235, 35, fill="#F3F6F9", stroke=LINE, sw=1.2, rx=16); text(1121, 758, "rank score • direction control", 14, GREY, "400", "middle", True)
tests=[("1", "Within-cohort coupling", "Spearman ρ • 2,000-bootstrap 95% CI", BLUE), ("2", "Across-cohort portability", "relationship matrix • study-level effects", TEAL), ("3", "Meta-analytic stability", "random-effects pooling • I² • leave-one-out", AMBER), ("4", "Competing explanations", "nested models • coefficients • adjusted R² • residual SD", PURPLE)]
for i,(n,t,s,col) in enumerate(tests):
    y=796+i*35; dot(600,y-5,col,17); text(600,y+2,n,15,"white","700","middle"); text(630,y-7,t,17,INK,"700"); text(630,y+13,s,13,GREY)
block(1120, 920, ["Unaccounted variation was not relabelled as", "RA-specific biology"], 12, GREY, "400", 17, "middle")

# Panel d
text(27, 989, "d", 28, INK, "700")
text(64, 989, "Prespecified interpretive rule and observed evidence pattern", 25, INK, "700")
rect(28, 1006, 1241, 147, fill="#F9FBFD", stroke=LINE, sw=2)
rect(50, 1025, 262, 105, fill=LIGHT_TEAL, stroke="#A8CDBF", sw=2, rx=13); text(181, 1070, "Measurable signal", 21, "#2C6B61", "700", "middle"); text(181, 1097, "continuous RA-derived", 15, INK, "400", "middle"); text(181, 1120, "projection in OA", 15, INK, "400", "middle")
line(315, 1078, 345, 1078, stroke=GREY, sw=4, marker="arrow-blue")
rect(356, 1025, 460, 120, fill="white", stroke=LINE, sw=2, rx=13); text(586, 1052, "What must travel with the score?", 18, INK, "700", "middle")
line(385, 1080, 740, 1080, AMBER, 5); [dot(x,1080,AMBER,10) for x in [423,521,628,715]]; text(563, 1103, "general inflammatory coupling — comparatively portable", 15, AMBER, "700", "middle")
line(385, 1120, 740, 1120, PURPLE, 5); [dot(x,1120,PURPLE,9) for x in [452,568,667,735]]; text(563, 1140, "APC/MHC-II • myeloid • IFN architecture — context dependent", 14, PURPLE, "700", "middle")
line(820, 1078, 850, 1078, stroke=GREY, sw=4, marker="arrow-blue")
rect(861, 1025, 386, 120, fill=LIGHT_AMBER, stroke="#D9BA67", sw=2, rx=13); text(1054, 1053, "Evidence boundary", 19, INK, "700", "middle"); text(1054, 1083, "Detectable resemblance", 18, RED, "700", "middle"); text(1054, 1106, "≠", 22, INK, "700", "middle"); text(1054, 1129, "stable RA-like OA subtype", 17, GREY, "700", "middle"); text(1054, 1145, "or shared disease identity", 12, GREY, "400", "middle")
text(31, 1194, "Molecular inference only: no clinical prediction, treatment response or causality was inferred. RA, rheumatoid arthritis; OA, osteoarthritis; APC, antigen-presenting-cell marker panel; IFN, interferon.", 12, GREY)
add('</svg>')

svg_path = OUT / "Figure1_IR_VECTOR_FINAL.svg"
svg_path.write_text("\n".join(svg), encoding="utf-8")
print(svg_path)
