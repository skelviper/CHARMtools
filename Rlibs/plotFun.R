pubReady <- function(p)
{
    # remove axis title for ggplot object
    library(ggthemes)
    p + theme_void(base_size = 15) + theme(title = element_blank())
}
