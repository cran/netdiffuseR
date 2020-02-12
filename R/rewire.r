#' Graph rewiring algorithms
#'
#' Changes the structure of a graph by altering ties.
#'
#' @inheritParams rgraph_ws
#' @templateVar undirected TRUE
#' @template graph_template
#' @param p Either a [0,1] vector with rewiring probabilities (\code{algorithm="endpoints"}),
#' or an integer vector with number of iterations (\code{algorithm="swap"}).
#' @param copy.first Logical scalar. When \code{TRUE} and \code{graph} is dynamic uses
#' the first slice as a baseline for the rest of slices (see details).
#' @param pr.change Numeric scalar. Probability ([0,1]) of doing a rewire (see details).
#' @param algorithm Character scalar. Either \code{"swap"}, \code{"endpoints"}, or \code{"qap"}
#' (see \code{\link{rewire_qap}}).
#' @param althexagons Logical scalar. When \code{TRUE} uses the compact alternating
#' hexagons algorithm (currently ignored [on development]).
#' @details
#' The algorithm \code{"qap"} is described in \code{\link{rewire_qap}}, and only
#' uses \code{graph} from the arguments (since it is simply relabelling the graph).
#'
#'
#' In the case of "swap" and "endpoints", both algorithms are implemented
#' sequentially, this is, edge-wise checking self edges and multiple edges over
#' the changing graph; in other words, at step
#' \eqn{m} (in which either a new endpoint or edge is chosen, depending on the algorithm),
#' the algorithms verify whether the proposed change creates either multiple edges
#' or self edges using the resulting graph at step \eqn{m-1}.
#'
#' The main difference between the two algorithms is that the \code{"swap"} algorithm
#' preserves the degree sequence of the graph and \code{"endpoints"} does not.
#' The \code{"swap"} algorithm is specially useful to asses the non-randomness of
#' a graph's structural properties, furthermore it is this algorithm the one used
#' in the \code{\link{struct_test}} routine implemented in \pkg{netdiffuseR}.
#'
#' Rewiring assumes a weighted network, hence \eqn{G(i,j) = k = G(i',j')},
#' where \eqn{i',j'} are the new end points of the edge and \eqn{k} may not be equal
#' to one.
#'
#' In the case of dynamic graphs, when \code{copy.first=TRUE}, after rewiring the
#' first slice--\eqn{t=1}--the rest of slices are generated by rewiring the rewired
#' version of the first slice. Formally:
#'
#' \deqn{%
#' G(t)' = \left\{\begin{array}{ll}
#' R(G(t)) & \mbox{if }t=1 \\
#' R(G(1)') & \mbox{otherwise}
#' \end{array}
#' \right.
#' }{%
#' G(t)' = R(G(t))  if t=1,
#'         R(G(1)') otherwise
#' }
#'
#' Where \eqn{G(t)} is the t-th slice, \eqn{G(t)'} is the t-th rewired slice, and
#' \eqn{R} is the rewiring function. Otherwise, \code{copy.first=FALSE} (default),
#' The rewiring function is simply \eqn{G(t)' = R(G(t))}.
#'
#' The following sections describe the way both algorithms were implemented.
#'
#' @section \emph{Swap} algorithm:
#' The \code{"swap"} algorithm chooses randomly two edges \eqn{(a,b)} and
#' \eqn{(c,d)} and swaps the 'right' endpoint of boths such that we get
#' \eqn{(a,d)} and \eqn{(c,b)} (considering self and multiple edges).
#'
#' Following Milo et al. (2004) testing procedure, the algorithm shows to be
#' well behaved in terms of been unbiased, so after each iteration each possible
#' structure of the graph has the same probability of been generated. The algorithm
#' has been implemented as follows:
#'
#' Let \eqn{E} be the set of edges of the graph \eqn{G}. For \eqn{i=1} to \eqn{p}, do:
#' \enumerate{
#'  \item With probability \code{1-pr.change} got to the last step.
#'  \item Choose \eqn{e0=(a, b)} from \eqn{E}. If \code{!self & a == b} then go to the last step.
#'  \item Choose \eqn{e1=(c, d)} from \eqn{E}. If \code{!self & c == d } then go to the last step.
#'  \item Define \eqn{e0'=(a, d)} and \eqn{e1' = (c, b)}. If \code{!multiple & [G[e0']!= 0 | G[e1'] != 0]} then go to the last step.(*)
#'  \item Define \eqn{v0 = G[e0]} and \eqn{v1 = G[e1]}, set \eqn{G[e0]=0} and \eqn{G[e1]=0}
#'  (and the same to the diagonally opposed coordinates in the case of undirected graphs)
#'  \item Set \eqn{G[e0'] = v0} and \eqn{G[e1'] = v1} (and so with the diagonally opposed coordinates
#'  in the case of undirected graphs).
#'  \item Next i.
#' }
#'
#' (*) When \code{althexagons=TRUE}, the algorithm changes and applies what Rao et al.
#' (1996) describe as Compact Alternating Hexagons. This modification assures the
#' algorithm to be able to achieve any structure. The algorithm consists on doing
#' the following swapping: \eqn{(i1i2,i1i3,i2i3,i2i1,i3i1,i3i2)} with values
#' \eqn{(1,0,1,0,1,0)} respectively with \eqn{i1!=i2!=i3}. See the examples and
#' references.
#'
#' In Milo et al. (2004) is suggested that in order for the rewired graph to be independent
#' from the original one researchers usually iterate around \code{nlinks(graph)*100}
#' times, so \code{p=nlinks(graph)*100}. On the other hand in Ray et al (2012)
#' it is shown that in order to achive such it is needed to perform
#' \code{nlinks(graph)*log(1/eps)}, where \code{eps}\eqn{\sim}1e-7, in other words,
#' around \code{nlinks(graph)*16}. We set the default to be 20.
#'
#' In the case of Markov chains, the variable \code{pr.change} allows making the
#' algorithm aperiodic. This is relevant only if the
#' probability self-loop to a particular state is null, for example, if
#' we set \code{self=TRUE} and \code{muliple=TRUE}, then in every step the
#' algorithm will be able to change the state. For more details see
#' Stanton and Pinar (2012) [p. 3.5:9].
#'
#'
#' @section \emph{Endpoints} algorithm:
#'
#' This reconnect either one or both of the endpoints of the edge randomly. As a big
#' difference with the swap algorithm is that this does not preserves the degree
#' sequence of the graph (at most the outgoing degree sequence). The algorithm is
#' implemented as follows:
#'
#' Let \eqn{G} be the baseline graph and \eqn{G'} be a copy of it. Then, For \eqn{l=1} to \eqn{|E|} do:
#'
#' \enumerate{
#'  \item Pick the \eqn{l}-th edge from \eqn{E}, define it as \eqn{e = (i,j)}.
#'  \item Draw \eqn{r} from \eqn{U(0,1)}, if \eqn{r > p} go to the last step.
#'  \item If \code{!undirected & i < j} go to the last step.
#'  \item Randomly select a vertex \eqn{j'} (and \eqn{i'} if \code{both_ends==TRUE}).
#'        And define \eqn{e'=(i, j')} (or \eqn{e'=(i', j')} if \code{both_ends==TRUE}).
#'  \item If \code{!self &} \code{i==j}' (or if \code{both_ends==TRUE & i'==j'}) go to the last step.
#'  \item If \code{!multiple & G'[e']!= 0} then go to the last step.
#'  \item Define \eqn{v = G[e]}, set \eqn{G'[e] = 0} and \eqn{G'[e'] = v} (and the
#'        same to the diagonally opposed coordinates in the case of undirected graphs).
#'  \item Next \eqn{l}.
#' }
#'
#' The endpoints algorithm is used by default in \code{\link{rdiffnet}} and used
#' to be the default in \code{\link{struct_test}} (now \code{swap} is the default).
#'
#' @references
#' Watts, D. J., & Strogatz, S. H. (1998). Collectivedynamics of "small-world" networks.
#' Nature, 393(6684), 440–442. \url{http://doi.org/10.1038/30918}
#'
#' Milo, R., Kashtan, N., Itzkovitz, S., Newman, M. E. J., & Alon, U.
#' (2004). On the uniform generation of random graphs with prescribed degree sequences.
#' Arxiv Preprint condmat0312028, cond-mat/0, 1–4. Retrieved from
#' \url{http://arxiv.org/abs/cond-mat/0312028}
#'
#' Ray, J., Pinar, A., and Seshadhri, C. (2012).
#' Are we there yet? When to stop a Markov chain while generating random graphs.
#' pages 1–21.
#'
#' Ray, J., Pinar, A., & Seshadhri, C. (2012). Are We There Yet? When to Stop a
#' Markov Chain while Generating Random Graphs. In A. Bonato & J. Janssen (Eds.),
#' Algorithms and Models for the Web Graph (Vol. 7323, pp. 153–164).
#' Berlin, Heidelberg: Springer Berlin Heidelberg.
#' \url{http://doi.org/10.1007/978-3-642-30541-2}
#'
#' A . Ramachandra Rao, R. J. and S. B. (1996). A Markov Chain Monte Carlo Method
#' for Generating Random ( 0 , 1 ) -Matrices with Given Marginals. The Indian
#' Journal of Statistics, 58, 225–242.
#'
#' Stanton, I., & Pinar, A. (2012). Constructing and sampling graphs with a
#' prescribed joint degree distribution. Journal of Experimental Algorithmics,
#' 17(1), 3.1. \url{http://doi.org/10.1145/2133803.2330086}
#'
#' @family simulation functions
#' @export
#' @author George G. Vega Yon
#' @examples
#' # Checking the consistency of the "swap" ------------------------------------
#'
#' # A graph with known structure (see Milo 2004)
#' n <- 5
#' x <- matrix(0, ncol=n, nrow=n)
#' x <- as(x, "dgCMatrix")
#' x[1,c(-1,-n)] <- 1
#' x[c(-1,-n),n] <- 1
#'
#' x
#'
#' # Simulations (increase the number for more precision)
#' set.seed(8612)
#' nsim <- 1e4
#' w <- sapply(seq_len(nsim), function(y) {
#'  # Creating the new graph
#'  g <- rewire_graph(x,p=nlinks(x)*100, algorithm = "swap")
#'
#'  # Categorizing (tag of the generated structure)
#'  paste0(as.vector(g), collapse="")
#' })
#'
#' # Counting
#' coded <- as.integer(as.factor(w))
#'
#' plot(table(coded)/nsim*100, type="p", ylab="Frequency %", xlab="Class of graph", pch=3,
#'  main="Distribution of classes generated by rewiring")
#'
#' # Marking the original structure
#' baseline <- paste0(as.vector(x), collapse="")
#' points(x=7,y=table(as.factor(w))[baseline]/nsim*100, pch=3, col="red")
#'
# ' # Compact Alternating Hexagons ----------------------------------------------
# ' x <- matrix(c(0,0,1,1,0,0,0,1,0), ncol=3, nrow=3)
# '
# ' set.seed(123)
# ' nsim <- 1e4
# ' w <- sapply(seq_len(nsim), function(y) {
# '  g <- rewire_graph(x,p=nlinks(x)*20, algorithm = "swap", althexagons=TRUE)
# '  paste0(as.vector(g), collapse="")
# ' })
# '
# ' # Counting
# ' coded <- as.integer(as.factor(w))
# '
# ' plot(table(coded)/nsim*100, type="p", ylab="Frequency %", xlab="Class of graph", pch=3,
# ' main="Distribution of classes generated by rewiring")
# '
# ' # Marking the original structure
# ' baseline <- paste0(as.vector(x), collapse="")
# ' points(x=7,y=table(as.factor(w))[baseline]/nsim*100, pch=3, col="red")
rewire_graph <- function(graph, p,
                         algorithm="endpoints",
                         both.ends=FALSE, self=FALSE, multiple=FALSE,
                         undirected=getOption("diffnet.undirected"),
                         pr.change= ifelse(self, 0.5, 1),
                         copy.first=TRUE, althexagons=FALSE) {

  # Checking undirected (if exists)
  checkingUndirected(graph)

  # althexagons is still on development
  if (althexagons) {
    althexagons <- FALSE
    warning("The option -althexagons- is still on development. So it has been set to FALSE.")
  }

  # Checking copy.first
  # if (missing(copy.first)) copy.first <- FALSE

  cls <- class(graph)
  out <- if ("dgCMatrix" %in% cls) {
    rewire_graph.dgCMatrix(graph, p, algorithm, both.ends, self, multiple, undirected, pr.change, althexagons)
  } else if ("list" %in% cls) {
    rewire_graph.list(graph, p, algorithm, both.ends, self, multiple, undirected, pr.change, copy.first, althexagons)
  } else if ("matrix" %in% cls) {
    rewire_graph.dgCMatrix(
      methods::as(graph, "dgCMatrix"), p, algorithm, both.ends, self, multiple, undirected, pr.change, althexagons)
  } else if ("diffnet" %in% cls) {
    rewire_graph.list(graph$graph, p, algorithm, both.ends, self, multiple,
                      graph$meta$undirected, pr.change, copy.first, althexagons)
  } else if ("array" %in% cls) {
    rewire_graph.array(graph, p, algorithm, both.ends, self, multiple, undirected, pr.change, copy.first, althexagons)
  } else stopifnot_graph(graph)

  # If diffnet, then it must return the same object but rewired, and change
  # the attribute of directed or not
  if (inherits(graph, "diffnet")) {
    graph$meta$undirected <- undirected
    graph$graph <- out
    return(graph)
  }

  attr(out, "undirected") <- FALSE

  return(out)
}

# @rdname rewire_graph
rewire_graph.list <- function(graph, p, algorithm, both.ends, self, multiple, undirected,
                              pr.change, copy.first, althexagons) {
  t   <- length(graph)
  out <- graph

  # Names
  tn <- names(graph)
  if (!length(tn)) tn <- 1:t
  names(out) <- tn

  # Checking p
  if (length(p)==1)
    p <- rep(p, t)

  for (i in 1:t) {

    # Copy replaces the first from 2 to T with 1
    j <- ifelse(copy.first, 1, i)

    out[[i]] <- if (algorithm == "endpoints")
      rewire_endpoints(out[[j]], p[i], both.ends, self, multiple, undirected)
    else if (algorithm == "swap")
      rewire_swap(out[[j]], p[i], self, multiple, undirected, pr.change) #, althexagons)
    else if (algorithm == "qap")
      rewire_qap(out[[j]])
    else stop("No such rewiring algorithm: ", algorithm)

    # Names
    rn <- rownames(graph[[i]])
    if (!length(rn)) rn <- 1:nrow(graph[[i]])
    dimnames(out[[i]]) <- list(rn, rn)
  }

  out
}

# @rdname rewire_graph
rewire_graph.dgCMatrix <- function(graph, p, algorithm, both.ends, self, multiple, undirected, pr.change, althexagons) {
  out <- if (algorithm == "endpoints")
    rewire_endpoints(graph, p, both.ends, self, multiple, undirected)
  else if (algorithm == "swap")
    rewire_swap(graph, p, self, multiple, undirected, pr.change) #, althexagons)
  else if (algorithm == "qap")
    rewire_qap(graph)
  else stop("No such rewiring algorithm: ", algorithm)

  rn <- rownames(out)
  if (!length(rn)) rn <- 1:nrow(out)
  dimnames(out) <- list(rn, rn)
  out
}

# @rdname rewire_graph
rewire_graph.array <-function(graph, p, algorithm, both.ends, self, multiple, undirected,
                              pr.change, copy.first, althexagons) {
  n   <- dim(graph)[1]
  t   <- dim(graph)[3]
  out <- apply(graph, 3, methods::as, Class="dgCMatrix")

  # Checking time names
  tn <- dimnames(graph)[[3]]
  if (!length(tn)) tn <- 1:t
  names(out) <- tn

  return(rewire_graph.list(out, p, algorithm, both.ends, self, multiple, undirected,
                    pr.change, copy.first, althexagons))
}

#' Permute the values of a matrix
#'
#' \code{permute_graph} Shuffles the values of a matrix either considering
#' \emph{loops} and \emph{multiple} links (which are processed as cell values
#' different than 1/0). \code{rewire_qap} generates a new graph \code{graph}\eqn{'}
#' that is isomorphic to \code{graph}.
#' @templateVar self TRUE
#' @templateVar multiple TRUE
#' @template graph_template
#' @author George G. Vega Yon
#' @return A permuted version of \code{graph}.
#' @examples
#' # Simple example ------------------------------------------------------------
#' set.seed(1231)
#' g <- rgraph_ba(t=9)
#' g
#'
#' # These preserve the density
#' permute_graph(g)
#' permute_graph(g)
#'
#' # These are isomorphic to g
#' rewire_qap(g)
#' rewire_qap(g)
#'
#' @references
#'
#' Anderson, B. S., Butts, C., & Carley, K. (1999). The interaction of size and
#' density with graph-level indices. Social Networks, 21(3), 239–267.
#' \url{http://dx.doi.org/10.1016/S0378-8733(99)00011-8}
#'
#' Mantel, N. (1967). The detection of disease clustering and a generalized
#' regression approach. Cancer Research, 27(2), 209–20.
#' \url{http://cancerres.aacrjournals.org/content/27/2_Part_1/209}
#'
#' @seealso This function can be used as null distribution in \code{struct_test}
#' @family simulation functions
#' @export
#' @aliases CUG QAP
permute_graph <- function(graph, self=FALSE, multiple=FALSE) {

  # Changing class
  cls <- class(graph)
  x <- if ("matrix" %in% cls) methods::as(graph, "dgCMatrix")
  else if ("list" %in% cls) lapply(graph, methods::as, Class="dgCMatrix")
  else if ("diffnet" %in% cls) graph$graph
  else if ("array" %in% cls) apply(graph, 3, methods::as, Class="dgCMatrix")
  else if ("dgCMatrix" %in% cls) graph
  else stopifnot_graph(graph)

  if (any(c("list", "array") %in% cls) & !("matrix" %in% cls)) {

    ans <- lapply(x, permute_graph_cpp, self=self, multiple=multiple)

  } else if ("diffnet" %in% cls) {
    ans <- graph
    ans$graph <- lapply(x, permute_graph_cpp, self=self, multiple=multiple)
  } else {
    ans <- permute_graph_cpp(x, self, multiple)
  }

  return(ans)

}

#' @export
#' @rdname permute_graph
rewire_permute <- permute_graph

#' @export
#' @rdname permute_graph
rewire_qap <- function(graph) {

  neword <- order(runif(nnodes(graph)))
  rewirefun <- function(graph) {
    graph[neword, neword]
  }

  # Changing class
  cls <- class(graph)
  x <- if ("matrix" %in% cls) methods::as(graph, "dgCMatrix")
  else if ("list" %in% cls) lapply(graph, methods::as, Class="dgCMatrix")
  else if ("diffnet" %in% cls) graph$graph
  else if ("array" %in% cls) apply(graph, 3, methods::as, Class="dgCMatrix")
  else if ("dgCMatrix" %in% cls) graph
  else
    stopifnot_graph(graph)

  if (any(c("diffnet", "list") %in% cls) | (("array" %in% cls) & length(dim(graph)) == 3L )) {

    ans <- lapply(x, rewirefun)

    if (inherits(graph, "diffnet")) {
      # Naming
      neword <- match(neword, nodes(graph))

      graph$graph <- ans
      graph$graph <- lapply(graph$graph, Matrix::unname)
      graph$meta$ids <- graph$meta$ids[neword]

      # Attributes
      if (nrow(graph$vertex.static.attrs)) {
        graph$vertex.static.attrs <- graph$vertex.static.attrs[neword,,drop=FALSE]
      }
      if (nrow(graph$vertex.dyn.attrs[[1]])) {
        graph$vertex.dyn.attrs <- lapply(graph$vertex.dyn.attrs, function(y) {
          y[neword,,drop=FALSE]
        })
      }

      # Adoptions
      graph$cumadopt <- graph$cumadopt[neword,,drop=FALSE]
      graph$adopt    <- graph$adopt[neword,,drop=FALSE]
      graph$toa      <- graph$toa[neword]

      return(graph)
    }


  } else {
    ans <- rewirefun(x)
  }

  return(ans)
}

