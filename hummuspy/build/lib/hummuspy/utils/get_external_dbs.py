import omnipath as om

def get_post_translational_data():
    ptm = om.interactions.PostTranslational()
    df = ptm.get()
    
    return df

