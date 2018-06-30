pgbart_predict <- function(x.test, model) {
  if(is.null(x.test))
    stop("ERROR: you miss a x.test parameter!")
  if(!is.matrix(x.test)){
    if(is.data.frame(x.test)){
      x.test <- as.matrix(makeind(x.test))
    }
    else{
      stop("ERROR: x.test must be a matrix or a data frame!")
    }
  }
  if(!is.character(model))
    stop("ERROR: model_file must be a character!")
  if(!file.exists(model))
    stop("ERROR: model_file does not exist!")

  metadata = as.numeric(read.table(model, nrows = 5)[[1]])
  rgy = as.double(metadata[1])
  ymin = as.double(metadata[2])
  ncpost = as.integer(metadata[3])
  binary = as.integer(metadata[4])
  binaryOffset = as.double(metadata[5])

  cres = .C('mpredict',
             as.integer(nrow(x.test)),
             as.integer(ncol(x.test)),
             as.double(x.test),
             as.character(model),
             tedraw = double(nrow(x.test)*ncpost))

  yhat.test = yhat.test.mean = NULL
  yhat.test = matrix(cres$tedraw, nrow = ncpost, byrow=T)
  if(!binary) {
    yhat.test = rgy * (yhat.test + .5) + ymin
    yhat.test.mean = apply(yhat.test, 2, mean)
  }
  if(binary) {
    yhat.test = yhat.test + binaryOffset
  }

  if(binary) {
    retval = list(
      call=match.call(),
      yhat.test=yhat.test,
      binaryOffset = binaryOffset
    )
  } else {
    retval = list(
      call=match.call(),
      yhat.test=yhat.test,
      yhat.test.mean=yhat.test.mean
    )
  }
  class(retval) = 'pgbart'
  return(invisible(retval))
}
