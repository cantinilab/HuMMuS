#!/bin/sh

INSTALL_DIR=$HOME/hummus_vm

#LOCAL_REPO="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
LOCAL_REPO="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; cd ../ >/dev/null 2>&1 ; pwd -P )"

LOCAL_IMAGE=$INSTALL_DIR/hummus.simg

LOCAL_IMAGE_SANDBOX=$INSTALL_DIR/sandbox

# display usage for current script
usage () {
    echo ""
    echo "Install the HUMMUS in $INSTALL_DIR."
    echo "Usage: $ bash $0 [-ueh] [-d path_to_CC_package]"
    echo ""
    echo "  -u      update hummus"
    echo "  -e      edit config file"
    echo "  -d      path to source code (development mode)"
    echo "  -h      print this help"
    echo ""
    exit 1
}

# compare versions
version_gt () {
    test "$(printf '%s\n' "$@" | sort -V | head -n 1)" != "$1";
}

# parse arguments
unset CREATE_IMAGE UPDATE_IMAGE EDIT_CONFIG CHANGE_ENV UPDATE_HUMMUS
CREATE_IMAGE=true
UPDATE_IMAGE=false
EDIT_CONFIG=false
CHANGE_ENV=false
UPDATE_HUMMUS=false
while getopts ':ued:hD' c
do
  case $c in
    u) CREATE_IMAGE=false; UPDATE_IMAGE=true; UPDATE_HUMMUS=true ;;
    e) CREATE_IMAGE=false; UPDATE_IMAGE=true; EDIT_CONFIG=true ;;
    d) CREATE_IMAGE=false; UPDATE_IMAGE=true; CHANGE_ENV=true; PATH_BRANCH=$OPTARG ;;
    D) DEBUG=true ;;
    h) usage ;; esac
done

# print variables if debugging
if [ "$DEBUG" = true ]
then
    echo NAME $NAME
    echo OPTARG $OPTARG
    echo CREATE_IMAGE $CREATE_IMAGE;
    echo UPDATE_IMAGE $UPDATE_IMAGE;
    echo EDIT_CONFIG $EDIT_CONFIG;
    echo CHANGE_ENV $CHANGE_ENV;
    echo UPDATE_HUMMUS $UPDATE_HUMMUS;
    echo PATH_BRANCH $PATH_BRANCH;
    exit 1;
fi

# check if valid path
if [ "$OPTARG" = "d" ] && [ ! -d "$PATH_BRANCH" ]
then
    printf -- "\033[31m ERROR: You need to specify a valid path when using the -d option. \033[0m\n";
    exit 1;
fi

# check if singularity is available
_=$(command -v singularity);
if [ "$?" != "0" ]
then
    printf -- "\033[31m ERROR: You don't seem to have Singularity installed \033[0m\n";
    printf -- 'Follow the guide at: https://www.sylabs.io/guides/2.6/user-guide/installation.html\n';
    exit 1;
fi;

# check singularity version
SINGULARITY_MIN_VERSION=2.5.0
SINGULARITY_INSTALLED_VERSION="$(singularity --version)"
if version_gt $SINGULARITY_MIN_VERSION $SINGULARITY_INSTALLED_VERSION
then
    printf -- "\033[31m ERROR: Update Singularity, we require at least version ${SINGULARITY_MIN_VERSION} (${SINGULARITY_INSTALLED_VERSION} detected) \033[0m\n";
    printf -- 'Follow the guide at: https://www.sylabs.io/guides/2.6/user-guide/installation.html\n';
    exit 2;
fi


# check if git is available
_=$(command -v git);
if [ "$?" != "0" ]
then
    printf -- "\033[31m ERROR: You don't seem to have Git installed \033[0m\n";
    printf -- 'Follow the guide at: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git\n';
    exit 3;
fi;




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
    sudo singularity build --sandbox $LOCAL_IMAGE_SANDBOX $SINGULARITY_DEFINITION  #cc_py37.def;
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

# check if a local singularity image is available, otherwise copy
if [ "$UPDATE_IMAGE" = true ]
then
    cd $INSTALL_DIR;
    # remove the previous image
    sudo rm -f $LOCAL_IMAGE;

    # change pythonpath in the image permanently to different CC repository
    # if [ "$CHANGE_ENV" = true ]
    # then
    #     text_replace_py="/export PYTHONPATH/c\\    export PYTHONPATH=\""$PATH_BRANCH"/package\":\$PYTHONPATH";
    #     text_replace_conf="/export CC_CONFIG/c\\    export CC_CONFIG=\""$PATH_BRANCH"/cc_config.json\"";
    #     sudo singularity exec  --writable $LOCAL_IMAGE_SANDBOX sed -i "$text_replace_py" /environment;
    #     sudo singularity exec  --writable $LOCAL_IMAGE_SANDBOX sed -i "$text_replace_conf" /environment;
    # fi

    # update CC to latest
    if [ "$UPDATE_HUMMUS" = true ]
    then
        # update sandbox
        sudo singularity exec  --writable $LOCAL_IMAGE_SANDBOX git --git-dir=/opt/chemical_checker/.git pull;
        if [ $? -eq 0 ]; then
            printf -- '\033[32m SUCCESS: Pulled latest Chemical Checker source code. \033[0m\n';
        else
            printf -- '\033[31m ERROR: Cannot update sandbox image. \033[0m\n';
            exit 7;
        fi
    fi

    # modify the config file
    if [ "$EDIT_CONFIG" = true ]
    then
        cd $INSTALL_DIR;
        sudo singularity exec  --writable $LOCAL_IMAGE_SANDBOX vi /opt/chemical_checker/setup/cc_config.json
        # generate image from sandbox
        sudo singularity build $LOCAL_IMAGE $LOCAL_IMAGE_SANDBOX
        if [ $? -eq 0 ]; then
            printf -- '\033[32m SUCCESS: Image file created. \033[0m\n';
        else
            printf -- '\033[31m ERROR: Cannot create image. \033[0m\n';
            exit 9;
        fi
    fi

    # generate image from sandbox
    sudo singularity build $LOCAL_IMAGE $LOCAL_IMAGE_SANDBOX;
    if [ $? -eq 0 ]; then
        printf -- '\033[32m SUCCESS: Image file created. \033[0m\n';
    else
        printf -- '\033[31m ERROR: Cannot create image. \033[0m\n';
        exit 8;
    fi
fi

