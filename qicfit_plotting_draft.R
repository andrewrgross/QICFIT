### Plotting QICFIT data -- Andrew R Gross
### This script loads in the coordinates of the nodes used in cytoscape, then outputs a draft figure.
### It can be used to test results during development.  Eventually, it will be refined to generate high quality graphics.

########################################################################
### Header
########################################################################

library(ggplot2)
library(RColorBrewer)
library(colorRamps)

########################################################################
### Palette
########################################################################

start = 50
end = 100
nColors = (end - start) * 2
my_palette <- rev(colorRampPalette(c("black","red","orange","yellow","white"))(n = nColors))
my_palette <- rev(heat.colors(nColors))

colors <- c("white","yellow","orange","red")

########################################################################
### Input
########################################################################

node.coord <- read.csv("Z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/templates for cytoscape/gtex-ggplot-full_node_coord.csv")
dummy.data <- read.csv("C:/Users/grossar/Desktop/dummy data.csv")

node.coord$value <- dummy.data$Average
node.coord$source <- as.character(node.coord$source)
node.coord$source[5] <- 'AN cingulate cortex'

df <- mtcars[, c("mpg", "cyl", "wt")]
df$cyl <- as.factor(df$cyl)

########################################################################
### Plot
########################################################################

ggplot(data = node.coord, aes(x = x.position, y = -y.position)) + 
  geom_point(aes(fill = value), size = 10, shape = 22) +
  scale_fill_gradientn(colors = c("white","white","yellow","orange", "red"), values = NULL, space = "Lab", na.value = "grey50", guide = "colourbar") +
  geom_text(aes(label = source),nudge_y = -30, size = 5) +
  geom_text(aes(label = value), size = 5.5, color = 'white') + 
  theme_minimal() +
  theme(panel.grid = element_blank(),axis.text = element_blank(),axis.title = element_blank(),legend.position = 'none')


theme(panel.grid.minor = element_line(colour = "red", linetype = "dotted"))


# Change colors and shapes manually
ggplot(df, aes(x=wt, y=mpg, group=cyl)) +
  geom_point(aes(shape=cyl, color=cyl), size=2)+
  scale_shape_manual(values=c(3, 16, 17))+
  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))+
  theme(legend.position="top")
# Change the point size manually
ggplot(df, aes(x=wt, y=mpg, group=cyl)) +
  geom_point(aes(shape=cyl, color=cyl, size=cyl))+
  scale_shape_manual(values=c(3, 16, 17))+
  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))+
  scale_size_manual(values=c(2,3,4))+
  theme(legend.position="top")