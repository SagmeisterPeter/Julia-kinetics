#using PlotlyJS
#using DataFrames
"""
    plot_paracoords

returns the paracoordinate plot of different input parameters 
"""

function plot_paracoords(df,string)
   mytrace = parcoords(line= attr(colorscale="telerose_r"),
        #;line = attr(color=["df.Ester"], colorscale=[(0,"red"), (0.5,"green"),(1,"blue")]),
        #color="Ester", colorscale="tealrose_r",
        dimensions = [ attr(range = [minimum(df.Temperature),maximum(df.Temperature)], label = "Temperature", values = df.Temperature), 
                       attr(range = [minimum(df.Total_Flow),maximum(df.Total_Flow)], label = "Total_Flow", values = df.Total_Flow),
                       attr(range = [minimum(df.Ester),maximum(df.Ester)], label = "Ester", values = df.Ester), 
                       attr(range = [minimum(df.Amine),maximum(df.Amine)], label = "Amine", values = df.Amine),
                       attr(range = [minimum(df.TBD),maximum(df.TBD)], label = "TBD", values = df.TBD),
                     #  attr(range = [minimum(tres),maximum(df.tres)], label = "tres", values = df.tres),
                       attr(range = [minimum(conversion),maximum(df.conversion)], label = "conversion", values = df.conversion),
                       attr(range = [minimum(productivity),maximum(df.productivity)], label = "productivity", values = df.productivity),
        ]);

    myplot = plot(mytrace)

    open("./"+string+".html", "w") do io
    PlotlyBase.to_html(io, myplot.plot)
   
    end

    return myplot
end


