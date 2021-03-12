#!/software/bin/perl -w

use strict;
use Cwd qw(cwd);

my $dir = cwd;

open my $prekt, '<', "sp2.kt"; 
my @pre = <$prekt>;
my @a = split "\t", $pre[0];
my $smallpre = $a[2];
my @b = split  "\t", $pre[-1];
my $bigpre = $b[2];

open my $postkt, '<', "sp1.kt"; 
my @post = <$postkt>;
my @c = split "\t", $post[0];
my $bigpost = $c[2];
my @d = split  "\t", $post[-1];
my $smallpost = $d[2];


open(OUT,">./circos.conf");

print OUT "\# circos.conf\n";
print OUT "\n";
print OUT "karyotype = ${dir}/sp1.kt,${dir}/sp2.kt\n";
print OUT "\n";
print OUT "\<ideogram\>\n";
print OUT "\n";
print OUT "\<spacing\>\n";
print OUT "default = 0.005r\n";
print OUT "\n";
print OUT "\# You can increase the spacing between specific ideograms.\n";
print OUT "\<pairwise ${bigpre} ${bigpost}\>\n";
print OUT "spacing = 20r\n";
print OUT "\</pairwise\>\n";
print OUT "\n";
print OUT "\<pairwise ${smallpost} ${smallpre}\>\n";
print OUT "spacing = 20r\n";
print OUT "\</pairwise\>\n";
print OUT "\n";
print OUT "\</spacing\>\n";
print OUT "\n";
print OUT "\# Ideogram position, thickness and fill.\n";
print OUT "radius           = 0.90r\n";
print OUT "thickness        = 20p\n";
print OUT "fill             = yes  \n";
print OUT "\#stroke_color     = dgrey\n";
print OUT "\#stroke_thickness = 2p   \n";
print OUT "\n";
print OUT "show_label     = yes\n";
print OUT "label_with_tag = yes\n";
print OUT "label_font     = light\n";
print OUT "label_radius   = dims(ideogram,radius_outer) + 0.1r\n";
print OUT "label_center   = yes\n";
print OUT "label_size     = 15p\n";
print OUT "label_color    = black\n";
print OUT "label_parallel = no\n";
print OUT "label_case     = upper \n";
print OUT "\#label_format   = eval(sprintf(\"chr%s\",var(label)))\n";
print OUT "\n";
print OUT "\</ideogram\>\n";
print OUT "\n";
print OUT "\n";
print OUT "\<links\>\n";
print OUT "\n";
print OUT "\<link\>\n";
print OUT "file =${dir}/links.txt\n";
print OUT "\n";
print OUT "radius = 0.99r\n";
print OUT "crest  = 1\n";
print OUT "ribbon           = yes\n";
print OUT "flat             = yes\n";
print OUT "\#stroke_color     = vdgrey\n";
print OUT "stroke_thickness = 1\n";
print OUT "\#color            = grey_a3\n";
print OUT "\n";
print OUT "bezier_radius        = 0r\n";
print OUT "bezier_radius_purity = 0.5\n";
print OUT "\n";
print OUT "\</link\>\n";
print OUT "\n";
print OUT "\</links>\n";
print OUT "\n";
print OUT "\################################################################\n";
print OUT "\# The remaining content is standard and required. It is imported from\n";
print OUT "\# default files in the Circos distribution.\n";
print OUT "\#\n";
print OUT "\# These should be present in every Circos configuration file and\n";
print OUT "\# overridden as required. To see the content of these files, \n";
print OUT "\# look in etc/ in the Circos distribution.\n";
print OUT "\#\n";
print OUT "\# It's best to include these files using relative paths. This way, the\n";
print OUT "\# files if not found under your current directory will be drawn from\n";
print OUT "\# the Circos distribution. \n";
print OUT "\#\n";
print OUT "\# As always, centralize all your inputs as much as possible.\n";
print OUT "\n";
print OUT "\<image\>\n";
print OUT "\# Included from Circos distribution.\n";
print OUT "\<\<include etc/image.conf\>\>\n";
print OUT "\</image\>\n";
print OUT "\n";
print OUT "\# RGB/HSV color definitions, color lists, location of fonts, fill\n";
print OUT "\# patterns. Included from Circos distribution.\n";
print OUT "\#\n";
print OUT "\# In older versions of Circos, colors, fonts and patterns were\n";
print OUT "\# included individually. Now, this is done from a central file. Make\n";
print OUT "\# sure that you're not importing these values twice by having\n";
print OUT "\#\n";
print OUT "\# *** DO NOT DO THIS ***\n";
print OUT "\# <colors>\n";
print OUT "\# <<include etc/colors.conf>>\n";
print OUT "\# <colors>\n";
print OUT "\# **********************\n";
print OUT "\<\<include etc/colors_fonts_patterns.conf\>\> \n";
print OUT "\n";
print OUT "\# Debugging, I/O an dother system parameters\n";
print OUT "\# Included from Circos distribution.\n";
print OUT "\<\<include etc/housekeeping.conf\>\>\n"; 


