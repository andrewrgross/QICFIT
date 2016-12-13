### Plotting QICFIT data -- Andrew R Gross


########################################################################
### Header
########################################################################

library(ggplot2)

########################################################################
### Functions
########################################################################


########################################################################
### Input
########################################################################

node.coord <- read.csv("Z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/templates for cytoscape/gtex-ggplot-full_node_coord.csv")

df <- mtcars[, c("mpg", "cyl", "wt")]
df$cyl <- as.factor(df$cyl)

########################################################################
### Plot
########################################################################

ggplot(data = node.coord, aes(x = x.position, y = -y.position)) + 
  geom_point(color = 'black', fill = 'white', size = 12, shape = 22) +
  geom_text(aes(label = source),nudge_y = -30, size = 3.2) +
  theme_minimal()


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