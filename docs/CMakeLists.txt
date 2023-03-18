# Documentation

find_package(Doxygen QUIET)

if(DOXYGEN_FOUND)
    find_program(SPHINX_BUILD sphinx-build)
    mark_as_advanced(FORCE SPHINX_BUILD)

    configure_file(Doxyfile.in Doxyfile)
    configure_file(conf.py.in conf.py)

    add_custom_target(doc DEPENDS ${PROJECT_BINARY_DIR}/docs/html/index.html)

    file(GLOB_RECURSE RST_FILES ${PROJECT_SOURCE_DIR}/docs/*.rst)
    set(N2_ERR_IMGS
        n2_err_c_p.png
        n2_err_c_v.png
        n2_err_h.png
        n2_err_k.png
        n2_err_mu.png
        n2_err_rho.png
        n2_err_s.png
        n2_err_u.png
        n2_err_w.png
    )

    add_custom_command(
        OUTPUT
            ${PROJECT_BINARY_DIR}/docs/html/index.html
        COMMAND
            ${SPHINX_BUILD} -c ${PROJECT_BINARY_DIR}/docs -b html -d ${PROJECT_BINARY_DIR}/docs/.doctrees ${PROJECT_SOURCE_DIR}/docs html
        DEPENDS
            ${PROJECT_BINARY_DIR}/docs/conf.py
            ${PROJECT_BINARY_DIR}/docs/xml/index.xml
            ${RST_FILES}
            ${N2_ERR_IMGS}
    )

    add_custom_command(
        OUTPUT
            ${PROJECT_BINARY_DIR}/docs/xml/index.xml
        COMMAND
            doxygen Doxyfile
        DEPENDS
            ${PROJECT_BINARY_DIR}/docs/Doxyfile
            ${PROJECT_SOURCE_DIR}/include/*.h
    )

    add_custom_command(
        OUTPUT
            ${N2_ERR_IMGS}
        COMMAND
            ${CMAKE_COMMAND} -E env "PYTHONPATH=${PROJECT_BINARY_DIR}/python/src:$ENV{PYTHONPATH}"
            python3 plot_err.py n2.yml
        WORKING_DIRECTORY
            ${PROJECT_SOURCE_DIR}/docs/pyplots
        DEPENDS
            ${PROJECT_SOURCE_DIR}/docs/pyplots/plot_err.py
            ${PROJECT_SOURCE_DIR}/docs/pyplots/n2.yml
    )

    add_custom_command(
        TARGET doc
        POST_BUILD
        COMMAND ;
        COMMENT
           "Open ${PROJECT_BINARY_DIR}/docs/html/index.html in your browser to view the documentation."
    )

endif()