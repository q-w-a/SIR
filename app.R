#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)

# Define UI for application that draws a histogram
ui <- fluidPage(
    withMathJax(),
    # Application title
    titlePanel("SIR Model"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            # sliderInput("bins",
            #             "Number of bins:",
            #             min = 1,
            #             max = 50,
            #             value = 30),
            sliderInput("i",
                        "Initial Infected:",
                        min = 1,
                        max = 1e5/2,
                        value = 2100),
            sliderInput("r",
                        "Initial Recovered:",
                        min = 1,
                        max =  1e5/2,
                        value = 2500),
            sliderInput("s",
                        "Initial Susceptible:",
                        min = 5000,
                        max = 1e5,
                        value = 95400),
            sliderInput("mu",
                        "$$\\mu, \\text{ the birth rate:}$$",
                        min = 0,
                        max = 1,
                        value = 0.09),
            sliderInput("beta",
                        "$$\\beta, \\text{ the transmission rate:}$$",
                        min = 0,
                        max = 1,
                        value = 0.2),
            sliderInput("gamma",
                        "$$\\gamma = \\frac{1}{d}, \\text{where } d \\text{ is duration of infectiousness:}$$",
                        min = 0,
                        max = 0.8,
                        value = 0.02),
            width = 4
        ),

        mainPanel(
            plotOutput("distPlot",height=700),
            plotOutput("formula",height=100),
            width =8

        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$formula <- renderPlot({

        mu <- .05
        beta <- 0.2
        gamma <- .02

        mu <- input$mu
        beta <- input$beta
        gamma <- input$gamma

        i <- input$i
        r <- input$r
        s <- input$s


        params = paste0("Parameters:\\;", "\\mu=", mu, ", \\beta = ", beta, ", \\gamma = ", gamma, "}$")
        text = paste0("$\\overset{Initial\\;Conditions:\\;",
                      "S_{initial} = ", s,
                      ", I_{initial}=", i,
                      ", R_{initial}=", r, "}{",
                      params)


        ggplot() +
            annotate("text", x = 4, y = 12, size=7, label = latex2exp::TeX(text)) +
            theme_void()

    }, height = 100)

    output$distPlot <- renderPlot({

        ####### initial conditions #########

        i <- input$i
        r <- input$r
        s <- input$s

        N <- i + r +s


        mu <- .05
        beta <- 0.2
        gamma <- .02

        mu <- input$mu
        beta <- input$beta
        gamma <- input$gamma


        s_rate <- function(S,I) { mu*N-(1*beta* S * I)/N - mu*S}
        i_rate <- function(S,I) {(beta * S * I)/N - (gamma+mu)*I}
        r_rate <- function(I, R) {gamma*I-mu*R}


        s_rates <- c()
        i_rates <- c()
        r_rates <- c()

        t_change <- 1


        r0 <- beta/(mu+gamma)

        eq <- tibble(S = ((gamma + mu) / beta )*N,
                     I = (mu/beta)*N*(beta/(gamma+mu) - 1),
                     R = N-S-I)


        time <- seq(1, 400, length = 1000)
        times <- c(0)
        for (t in time) {
            ind <- length(s)
            times <- c(times, t)


            s_new <- s[ind] + s_rate(s[ind], i[ind]) * t_change
            i_new <- i[ind] + i_rate(s[ind], i[ind]) * t_change
            r_new <- r[ind] + r_rate(i[ind],r[ind]) * t_change

            s <- c(s, s_new)
            i <- c(i, i_new)
            r <- c(r, r_new)


            if (r0 <=1) { if(abs(s_new - s[ind] <= 1)) break}
            if (r0 > 1) { if (abs(s_new - eq$S) <= 1) break}

        }


        results <- tibble::tibble(time =times,
                                  susceptible = s,
                                  infected = i,
                                  recovered = r
        )

        custom_title <- paste0("$R_0 = \\frac{\\beta}{\\mu + \\gamma} =$", round(beta/(mu+gamma),3))

        s_star <- ((gamma + mu)/beta) * N
        i_star <- (mu/beta)*N*((beta)/(gamma+mu) - 1)
        r_star <- N - s_star - i_star


        end_eq <- paste0("Endemic Equilibrium at $(S^*,I^*,R^*) = \\frac{\\gamma + \\mu}{\\beta} N,",
                         "\\; \\frac{\\mu}{\\beta} N ( \\frac{\\beta}{\\gamma+\\mu} - 1), \\; N - S^* - I^*)$")
        end_eq_pt <- paste0("$ = (", round(eq$S,0),  ", ", round(eq$I,0),  ", ", round(eq$R,0), ")$")

        endemic <- ifelse(beta/(mu+gamma) > 1, paste0(end_eq, end_eq_pt), "")


        p <- results %>%
            tidyr::pivot_longer(c(susceptible,
                                  infected,
                                  recovered),
                                names_to = "category") %>%
            mutate(eq = case_when(
                category == "infected" ~ eq$I,
                category == "susceptible"~ eq$S,
                category == "recovered"~eq$R)) %>%
            ggplot(aes(x = time,
                       y = value,
                       color = category)) +
            labs(color = "",
                 title = custom_title) +
            geom_line(size = 1.3) +
            theme_bw() +
            theme(axis.title = element_text(size = 20, face = "bold"),
                  axis.text = element_text(size = 15),
                  legend.text = element_text(size = 20),
                  plot.title = element_text(size = 30, hjust = .5),
                  plot.caption = element_text(size = 13)) +
            scale_y_continuous(labels = scales::comma) +
            labs(x= "Time", y = "Number of People",
                 title = latex2exp::TeX(custom_title),
                 caption = latex2exp::TeX(endemic)) +
            scale_x_continuous(n.breaks = 10) +
            guides(color = guide_legend(override.aes = list(size = 4))) +
            scale_color_manual(values = c("#42AFB2", "#167134", "#BFA550"))

        if(r0<=1) { print(p)
        } else {
            p <- p +
            geom_hline(aes(yintercept = eq,
                           color = category),
                       linetype = 2)
        print(p)}
    }, height = 650)
}

# Run the application
shinyApp(ui = ui, server = server)
