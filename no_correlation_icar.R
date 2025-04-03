ws10 <- read.csv("pre_icar_data/y_ws10.csv", header=FALSE)
tp <- read.csv("pre_icar_data/y_tp.csv", header = FALSE)

node1 <- read.csv("pre_icar_data/node1.csv", header = FALSE)
node2 <- read.csv("pre_icar_data/node2.csv", header = FALSE)

node1 <- as.integer(node1[[1]])  # Convert node1 to a vector
node2 <- as.integer(node2[[1]])  # Convert node2 to a vector
y1 <- as.vector(ws10[[1]])  # Convert ws10 to a numeric vector
y2 <- as.vector(tp[[1]])  # Convert ws10 to a numeric vector