#!/bin/sh

INSTALL_DIR=$HOME/hummus_vm

#LOCAL_REPO="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
LOCAL_REPO="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; cd ../ >/dev/null 2>&1 ; pwd -P )"

LOCAL_IMAGE=$INSTALL_DIR/hummus.simg

LOCAL_IMAGE_SANDBOX=$INSTALL_DIR/sandbox

# parse arguments
unset CREATE_IMAGE 
CREATE_IMAGE=true

# if we are creating the image we remove image and sandbox
if [ "$CREATE_IMAGE" = true ]
then
    if [ -d "$INSTALL_DIR" ]
    then
        printf -- '\033[33m WARNING: %s already exists, delete it to proceed. \033[0m\n' $INSTALL_DIR;
        exit 0;
    fi
    mkdir $INSTALL_DIR;
    cd $INSTALL_DIR;
    SINGULARITY_DEFINITION="$LOCAL_REPO/container/recipe_singularity.def"

    RDEPSCRIPT="$LOCAL_REPO/container/install_rdeps.R"
    ENVCONFIG="$LOCAL_REPO/container/environment.yml"
    cp $SINGULARITY_DEFINITION $INSTALL_DIR;  # Node2vec binary
    mkdir -p ./container/;
    cp -r $RDEPSCRIPT ./container/;
    cp -r $ENVCONFIG ./container/;
    #cp -r $SNAP_MAKEFILE ./container/snap_makefile.config;

    printf -- 'Removing old singularity image...\n';
    sudo rm -f $LOCAL_IMAGE;
    sudo rm -rf $LOCAL_IMAGE_SANDBOX;

    printf -- 'Creating singularity sandbox image... \n';
    sudo singularity build --sandbox $LOCAL_IMAGE_SANDBOX $SINGULARITY_DEFINITION
    if [ $? -eq 0 ]; then
        printf -- '\033[32m SUCCESS: Image sandbox created correctly. \033[0m\n';
    else
        printf -- '\033[31m ERROR: Cannot create sandbox image. \033[0m\n';
        exit 5;
    fi

    # generate image from sandbox
    sudo singularity build $LOCAL_IMAGE $LOCAL_IMAGE_SANDBOX;
    if [ $? -eq 0 ]; then
        printf -- '\033[32m SUCCESS: Image file created. \033[0m\n';
    else
        printf -- '\033[31m ERROR: Cannot create image. \033[0m\n';
        exit 6;
    fi


fi

