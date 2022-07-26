def color_mapper(name):
    
    if 'katie' in name and 'ENCODE' in name:
        color = 'rosybrown'
    elif 'NATURE' in name:
        color = 'tomato'
    elif 'katie' in name:
        color = 'royalblue'
        if 'V5' in name:
            color = 'deepskyblue'
    elif  '262' in name.upper() or 'ENCODE' in name.upper():
        print(name.upper())
        if 'down' in name or 'ENCODE' in name:
            color = 'hotpink'
        else:
            color = 'pink'
    else:
        color = 'gold'
    return color