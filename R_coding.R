####### deal with the error in R
# tryCatch has a slightly complex syntax structure. However, once we understand the 4 parts which constitute a complete tryCatch call as shown below, it becomes easy to remember:

# expr: [Required] R code(s) to be evaluated
# error : [Optional] What should run if an error occured while evaluating the codes in expr
# warning : [Optional] What should run if a warning occured while evaluating the codes in expr
# finally : [Optional] What should run just before quitting the tryCatch call, irrespective of if expr ran successfully, with an error, or with a warning

tryCatch(
    expr = {
        # Your code...
        # goes here...
        # ...
    },
    error = function(e){ 
        # (Optional)
        # Do this if an error is caught...
    },
    warning = function(w){
        # (Optional)
        # Do this if an warning is caught...
    },
    finally = {
        # (Optional)
        # Do this at the end before quitting the tryCatch structure...
    }
)
