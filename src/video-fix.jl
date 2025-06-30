# add manually in Jupyter notebooks
#__precompile__(false)
#=
function Base.show(io::IO, ::MIME"text/html", vs::VideoStream)
    mktempdir() do dir
        path = save(joinpath(dir, "video.mp4"), vs)
        # <video> only supports infinite looping, so we loop forever even when a finite number is requested
        loopoption = "" # vs.options.loop ≥ 0 ? "loop" : ""
        Base.display(
            HTML("""<video style="max-width: 100%;" autoplay controls $loopoption><source src="data:video/x-m4v;base64,""" *
            Makie.Base64.base64encode(open(read, path)) *
            """" type="video/mp4"></video>""")
        )
    end
end
=#