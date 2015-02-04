import sys

from weblogolib import *
from weblogolib.colorscheme import *

# read argument
filename = sys.argv[1]

RAA = ColorScheme([
  ColorGroup( "DE",  "green",   "gp1"),
  ColorGroup( "NQ",      "purple", "gp2"),
  ColorGroup( "FWY",     "blue",   "gp3"),
  ColorGroup( "RHK",      "red",    "gp4"),
  ColorGroup("AILMV", "black",  "gp5"),
  ColorGroup("GP", "orange",  "gp6"),
  ColorGroup("ST", "yellow",  "gp7"),
  ColorGroup("C", "cyan",  "gp8"), ]
  )
test = ColorScheme([
  ColorGroup( "A",  "green",   "gp1"),
  ColorGroup( "B",      "purple", "gp2"),
  ColorGroup( "C",     "blue",   "gp3"),
  ColorGroup( "D",      "red",    "gp4"),
  ColorGroup("E", "black",  "gp5"),
  ColorGroup("F", "orange",  "gp6"),
  ColorGroup("G", "yellow",  "gp7"),
  ColorGroup("H", "cyan",  "gp8"), ]
  )
fin = open("data/" + filename + ".txt", "r")
seqs = read_seq_data(fin, alphabet="DENQFWYRHKAILMVGPSTC")
fin.close()
data = LogoData.from_seqs(seqs)
options = LogoOptions()
options.color_scheme = RAA
options.stack_width=40
format = LogoFormat(data, options)
eps = pdf_formatter( data, format)
fout = open("plots/" + filename + ".pdf", "wb")
fout.write(eps)
fout.close()
