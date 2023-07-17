library(igraph)

g <- make_empty_graph()
g <- make_graph(edges = c(1,2, 1,5), n=10, directed = FALSE)

plot(g)


g <- make_graph('Zachary')
plot(g)

# add mais vertices
g <- add_vertices(g, 3)
plot(g)

# add mais arestas
g <- add_edges(g, edges = c(1,35, 1,36, 34,37))
plot(g)

# operador plus +
g <- g + edges(1,35, 1,36, 34,37)
plot(g)

# encademaneto
g <- g %>% add_edges(edges=c(1,34)) %>% add_vertices(3) %>%
  add_edges(edges=c(38,39, 39,40, 40,38, 40,37))
g
plot(g)

get.edge.ids(g, c(1,34))
g <- delete_edges(g, 82)

g <- make_ring(10) %>% delete_edges("10|1")
plot(g)

g <- make_ring(5)
g <- delete_edges(g, get.edge.ids(g, c(1,5, 4,5)))
plot(g)

# grafo cordal
g1 <- graph_from_literal(A-B:C:I, B-A:C:D, C-A:B:E:H, D-B:E:F,
                         E-C:D:F:H, F-D:E:G, G-F:H, H-C:E:G:I,
                         I-A:H)
plot(g1)

is_chordal(g1, fillin=TRUE)

chordal_graph <- add_edges(g1, is_chordal(g1, fillin=TRUE)$fillin)
plot(chordal_graph)

# -- construindo gráficos
## mesmo grafico
graph1 <- make_tree(127, 2, mode = "undirected")
graph2 <- make_tree(127, 2, mode = "undirected")
identical_graphs(graph1,graph2)
plot(graph1)

## graficos diferentes
graph1 <- sample_grg(100, 0.2)
graph2 <- sample_grg(100, 0.2)
identical_graphs(graph1, graph2)

isomorphic(graph1, graph2)


## configurando e recuperando atributos
g <- make_graph(~ Alice-Bob:Claire:Frank, Claire-Alice:Dennis:Frank:Esther,
                George-Dennis:Frank, Dennis-Esther)
plot(g)

V(g)$age <- c(25, 31, 18, 23, 47, 22, 50) 
V(g)$gender <- c("f", "m", "f", "m", "m", "f", "m")
E(g)$is_formal <- c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE)
summary(g)

### outra forma
g <- make_graph(~ Alice-Bob:Claire:Frank, Claire-Alice:Dennis:Frank:Esther,
                George-Dennis:Frank, Dennis-Esther) %>%
  set_vertex_attr("age", value = c(25, 31, 18, 23, 47, 22, 50)) %>%
  set_vertex_attr("gender", value = c("f", "m", "f", "m", "m", "f", "m")) %>%
  set_edge_attr("is_formal", value = c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE))
summary(g)

### atribuir ou modificar um atributo
E(g)$is_formal

E(g)$is_formal[1] <- TRUE
E(g)$is_formal

g$date <- c("2023-07-17")
graph_attr(g, "date")

match(c("George"), V(g)$name)

V(g)$name[1:3] <- c("Alejandra", "Bruno", "Carmina")
V(g)

### deletar atributos
g <- delete_vertex_attr(g, "gender")
V(g)$gender

## propriedades estruturais de graficos
degree(g)
degree(g, 7)
degree(g, v=c(3,4,5))
degree(g, v=c("Carmina", "Frank", "Dennis"))
degree(g, "Bruno")

### intermediacao de borda
edge_betweenness(g)
ebs <- edge_betweenness(g)
as_edgelist(g)[ebs == max(ebs), ]

## consultando vertices e arestas com base em atributos
plot(g)
which.max(degree(g))

### ids impares
graph <- graph.full(n=10)
only_odd_vertices <- which(V(graph)%%2==1)
length(only_odd_vertices)

#### nomes dos indivíduos com menos de 30 anos
V(g)[age < 30]$name

### selecionando arestas
E(g)[.from(3)] # originarias de Carmina
E(g)[.from("Carmina")]

V(g)$gender <- c("f", "m", "f", "m", "m", "f", "m")
men <- V(g)[gender == "m"]$name
men

women <- V(g)[gender == "f"]$name
women

E(g)[men %--% women]

## tratando um grafo como uma matriz de adjacencia
get.adjacency(g)

## layouts e plotagem
layout <- layout_with_kk(g)
layout <- layout_as_tree(g, root = 2)

## desenhando um grafico usando um layout
layout <- layout_with_kk(g)
plot(g, layout = layout, main = "Social network with the Kamada-Kawai layout algorithm")

plot(g, layout = layout_with_fr,
     main = "Social network with the Fruchterman-Reingold layout algorithm")

V(g)$color <- ifelse(V(g)$gender == "m", "yellow", "red")
plot(g, layout = layout, vertex.label.dist = 3.5,
     main = "Social network - with genders as colors")

plot(g, layout=layout, vertex.label.dist=3.5, vertex.color=as.factor(V(g)$gender))

plot(g, layout=layout, vertex.label.dist=3.5, vertex.size=20,
     vertex.color=ifelse(V(g)$gender == "m", "yellow", "red"),
     edge.width=ifelse(E(g)$is_formal, 5, 1))
