library(heatmaply)
library(htmlwidgets)

input = commandArgs(trailingOnly=TRUE)

gridf = input[1]
#setwd('/home/alex/git/PlasCliques/')
#gridf = "claire_all.grid"

#cats = input[4:length(input)]
cats = c("Unannotated", "Rep", "Mob", "Partitioning", "Conjugation")
ncats = length(cats)
colourz = viridis(ncats)

bbs = read.table(file=gridf, header = FALSE, row.names = 1, sep = "\t")

CONST = 0
# Try const = 0, if more than X cols from dim then up until <X?
while (dim(bbs)[2]>400) {
  bbs = bbs[, colSums(bbs > 0) > CONST]
  row.order = hclust( dist(bbs), method = "ward.D" )$order
  bbs = bbs[row.order,]
  # Columns must have more than CONST values in a row, (needs cluster/initial row reordering)
  include = apply(bbs,2,function(x) max(rle(x>0)$lengths[rle(x>0)$values])>CONST)
  bbs = bbs[,include]
  CONST = CONST+1
}

print(paste0('CONST: ',CONST-1))


# If range 0-1 no legend, else get legend titles as input (rep,relaxase etc)

hm = heatmaply(bbs, showticklabels = c(FALSE,TRUE), hide_colorbar = TRUE,
               dendrogram = 'both', dist_method = 'binary', hclust_method = 'ward.D',
               label_names = c("<b>Isolate</b>","<b>Clique</b>","<i>Value</i>"),
               colors = c('white',colourz), grid_gap = 0.2, fontsize_row = 10, height = 600, plot_method = "ggplot")

js = "document.getElementsByClassName('plotly')[0].on('plotly_doubleclick', function(){Plotly.relayout(document.getElementsByClassName('plotly')[0], {'xaxis.autorange': true,'yaxis.autorange': true});document.querySelector('[data-title=\"Reset axes\"]').click();});"
p3 = prependContent(hm, onStaticRenderComplete(js), data = list(''))

# If cats input, legend...
legfn = tempfile()
svg(filename=legfn,bg=NA)
plot(100,100, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("center", cats, horiz = FALSE, fill = colourz, cex = 1)
dev.off()
leg = base64enc::base64encode(legfn)
tag = paste0("data:image/svg+xml;base64,",leg)
js = paste0("var node = document.createElement('style'); node.innerHTML = '#overlay { pointer-events: none; width: 100%; height: 100%; top: 0; left: 0; right: 0; bottom: 0; background-image: url(\\'",tag,"\\'); position: fixed; background-repeat: no-repeat; background-position: bottom -250px right -210px;}'; document.body.appendChild(node); zz = document.createElement('div'); zz.setAttribute('id', 'overlay'); document.body.appendChild(zz);")
p3 = prependContent(p3, onStaticRenderComplete(js), data = list(''))

saveWidget(p3, selfcontained = FALSE, file=paste0(gridf,'.html'))
