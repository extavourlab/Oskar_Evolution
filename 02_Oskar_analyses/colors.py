import seaborn as sns

class colors:
    def __init__(self):
        self.mapping_colors_r = {
                                'Diptera':"#fffbc8",
                                'Mecoptera':"#ffedb4",
                                'Lepidoptera':"#f9cdac",
                                'Trichoptera':"#f3aca2",
                                'Coleoptera':"#ee8b97",
                                'Neuroptera':"#dd889a",
                                'Hymenoptera':"#e96a8d",
                                'Psocoptera':"#db5087",
                                'Thysanoptera':"#b8428c",
                                'Blattodea':"#973490",
                                'Phasmatodea':"#742796",
                                'Orthoptera':"#5e1f88",
                                'Plecoptera':"#4d1a70",
                                'Ephemeroptera':"#3d1459",
                                'Zygentoma':"#2d0f41"
                                }
        
        self.legend_order = ['Zygentoma',
                    'Ephemeroptera',
                    'Plecoptera',
                    'Orthoptera',
                    'Phasmatodea',
                    'Blattodea',
                    'Thysanoptera',
                    'Psocoptera',
                    'Hymenoptera',
                    'Neuroptera', 
                    'Coleoptera',
                    'Trichoptera',
                    'Lepidoptera',
                    'Mecoptera',
                    'Diptera']
        self.light_gray = "#cccccc"
        self.gray = "#999999"
        self.dark_gray = "#555555"
        
    def order(self, name):
        if name in self.mapping_colors_r:
            return self.mapping_colors_r[name]
        else:
            return self.light_gray
        
    def get_pallette(self):
        res = []
        for order in self.legend_order:
            res.append(self.order(order))
        return sns.color_palette(res)