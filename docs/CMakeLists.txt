# Documentation

find_package(Doxygen QUIET)

if(DOXYGEN_FOUND)
    # List of all images with error evaluated against coolprop
    set(ERR_IMAGES "")

    macro(add_fluid NAME YAML_FILE)
        set(${NAME}_ERR_IMGS
            ${NAME}_err_c_p.png
            ${NAME}_err_c_v.png
            ${NAME}_err_h.png
            ${NAME}_err_k.png
            ${NAME}_err_mu.png
            ${NAME}_err_rho.png
            ${NAME}_err_s.png
            ${NAME}_err_u.png
            ${NAME}_err_w.png
        )

        add_custom_command(
            OUTPUT
                ${${NAME}_ERR_IMGS}
            COMMAND
                ${CMAKE_COMMAND} -E env "PYTHONPATH=${PROJECT_BINARY_DIR}/python/src:$ENV{PYTHONPATH}"
                python3 plot_err.py ${YAML_FILE}
            WORKING_DIRECTORY
                ${PROJECT_SOURCE_DIR}/docs/pyplots
            DEPENDS
                ${PROJECT_SOURCE_DIR}/docs/pyplots/plot_err.py
                ${PROJECT_SOURCE_DIR}/docs/pyplots/${YAML_FILE}
        )

        list(APPEND ERR_IMAGES ${${NAME}_ERR_IMGS})
    endmacro()

    find_program(SPHINX_BUILD sphinx-build)
    mark_as_advanced(FORCE SPHINX_BUILD)

    configure_file(Doxyfile.in Doxyfile)
    configure_file(conf.py.in conf.py)

    add_custom_target(doc DEPENDS ${PROJECT_BINARY_DIR}/docs/html/index.html)

    file(GLOB_RECURSE RST_FILES ${PROJECT_SOURCE_DIR}/docs/*.rst)
    add_fluid(air air.yml)
    add_fluid(he he.yml)
    add_fluid(n2 n2.yml)
    add_fluid(co2 co2.yml)
    add_fluid(o2 o2.yml)
    add_fluid(ch4 ch4.yml)

    add_custom_command(
        OUTPUT
            ${PROJECT_BINARY_DIR}/docs/html/index.html
        COMMAND
            ${SPHINX_BUILD} -c ${PROJECT_BINARY_DIR}/docs -b html -d ${PROJECT_BINARY_DIR}/docs/.doctrees ${PROJECT_SOURCE_DIR}/docs html
        DEPENDS
            ${PROJECT_BINARY_DIR}/docs/conf.py
            ${PROJECT_BINARY_DIR}/docs/xml/index.xml
            ${RST_FILES}
            ${ERR_IMAGES}
    )

    add_custom_command(
        OUTPUT
            ${PROJECT_BINARY_DIR}/docs/xml/index.xml
        COMMAND
            doxygen Doxyfile
        DEPENDS
            ${PROJECT_BINARY_DIR}/docs/Doxyfile
            ${PROJECT_SOURCE_DIR}/include/fprops/*.h
    )

    add_custom_command(
        TARGET doc
        POST_BUILD
        COMMAND ;
        COMMENT
           "Open ${PROJECT_BINARY_DIR}/docs/html/index.html in your browser to view the documentation."
    )

endif()
