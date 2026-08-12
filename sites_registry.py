"""Master Database of Training Sites for Multi-Site AI Training.

This module provides a registry of 40 historic and reference sites across Kerala/South India:
- 10 Landslide disaster sites (Sentinel-2 NDVI based label fetching)
- 10 Flood disaster sites (Sentinel-1 SAR based label fetching)
- 20 Safe baseline sites (Flat, urban, or geologically stable control sites)
"""

from typing import Dict, Any, Tuple

TRAINING_SITES: Dict[str, Dict[str, Any]] = {
    # =========================================================================
    # 1. LANDSLIDE DISASTER SITES (10)
    # =========================================================================
    'meppadi': {
        'bbox': (76.10, 11.50, 76.22, 11.60),
        'category': 'landslide',
        'label_type': 'sentinel2_ndvi',
        'pre_event': ('2024-01-01', '2024-06-30'),
        'post_event': ('2024-08-01', '2024-12-31'),
        'description': '2024 Wayanad landslide (Chooralmala/Mundakkai)',
    },
    'kavalappara': {
        'bbox': (76.10, 11.10, 76.20, 11.20),
        'category': 'landslide',
        'label_type': 'sentinel2_ndvi',
        'pre_event': ('2019-01-01', '2019-07-15'),
        'post_event': ('2019-09-01', '2019-12-31'),
        'description': '2019 Malappuram landslide (Kavalappara)',
    },
    'pettimudi': {
        'bbox': (77.00, 10.03, 77.10, 10.13),
        'category': 'landslide',
        'label_type': 'sentinel2_ndvi',
        'pre_event': ('2020-01-01', '2020-07-15'),
        'post_event': ('2020-09-01', '2020-12-31'),
        'description': '2020 Idukki landslide (Pettimudi)',
    },
    'koottickal': {
        'bbox': (76.80, 9.53, 76.90, 9.63),
        'category': 'landslide',
        'label_type': 'sentinel2_ndvi',
        'pre_event': ('2021-01-01', '2021-09-30'),
        'post_event': ('2021-11-01', '2022-03-31'),
        'description': '2021 Kottayam landslide (Koottickal)',
    },
    'puthumala': {
        'bbox': (76.17, 11.50, 76.27, 11.60),
        'category': 'landslide',
        'label_type': 'sentinel2_ndvi',
        'pre_event': ('2019-01-01', '2019-07-15'),
        'post_event': ('2019-09-01', '2019-12-31'),
        'description': '2019 Wayanad landslide (Puthumala)',
    },
    'rajamala': {
        'bbox': (77.01, 10.07, 77.11, 10.17),
        'category': 'landslide',
        'label_type': 'sentinel2_ndvi',
        'pre_event': ('2020-01-01', '2020-07-15'),
        'post_event': ('2020-09-01', '2020-12-31'),
        'description': '2020 Idukki landslide (Rajamala)',
    },
    'adimali': {
        'bbox': (76.91, 9.97, 77.01, 10.07),
        'category': 'landslide',
        'label_type': 'sentinel2_ndvi',
        'pre_event': ('2021-01-01', '2021-09-30'),
        'post_event': ('2021-11-01', '2022-03-31'),
        'description': '2021 Idukki debris flow (Adimali)',
    },
    'kokkayar': {
        'bbox': (76.90, 9.55, 77.00, 9.65),
        'category': 'landslide',
        'label_type': 'sentinel2_ndvi',
        'pre_event': ('2021-01-01', '2021-09-30'),
        'post_event': ('2021-11-01', '2022-03-31'),
        'description': '2021 Idukki landslide (Kokkayar)',
    },
    'vilangad': {
        'bbox': (75.80, 11.40, 75.90, 11.50),
        'category': 'landslide',
        'label_type': 'sentinel2_ndvi',
        'pre_event': ('2019-01-01', '2019-07-15'),
        'post_event': ('2019-09-01', '2019-12-31'),
        'description': '2019 Kozhikode landslide (Vilangad)',
    },
    'upputhode': {
        'bbox': (76.90, 9.80, 77.00, 9.90),
        'category': 'landslide',
        'label_type': 'sentinel2_ndvi',
        'pre_event': ('2022-01-01', '2022-07-15'),
        'post_event': ('2022-09-01', '2022-12-31'),
        'description': '2022 Idukki landslide (Upputhode)',
    },

    # =========================================================================
    # 2. FLOOD DISASTER SITES (10)
    # =========================================================================
    'chellanam': {
        'bbox': (76.28, 9.78, 76.35, 9.83),
        'category': 'flood',
        'label_type': 'sentinel1_sar',
        'pre_event': ('2018-07-01', '2018-07-25'),
        'post_event': ('2018-08-15', '2018-08-25'),
        'description': '2018 Chellanam coastal flood inundation',
    },
    'aluva': {
        'bbox': (76.31, 10.06, 76.41, 10.16),
        'category': 'flood',
        'label_type': 'sentinel1_sar',
        'pre_event': ('2018-07-01', '2018-07-25'),
        'post_event': ('2018-08-15', '2018-08-25'),
        'description': '2018 Aluva flood (Periyar river basin)',
    },
    'chengannur': {
        'bbox': (76.57, 9.27, 76.67, 9.37),
        'category': 'flood',
        'label_type': 'sentinel1_sar',
        'pre_event': ('2018-07-01', '2018-07-25'),
        'post_event': ('2018-08-15', '2018-08-25'),
        'description': '2018 Chengannur flood (Pampa river basin)',
    },
    'kuttanad': {
        'bbox': (76.38, 9.37, 76.48, 9.47),
        'category': 'flood',
        'label_type': 'sentinel1_sar',
        'pre_event': ('2018-07-01', '2018-07-25'),
        'post_event': ('2018-08-15', '2018-08-25'),
        'description': '2018 Kuttanad backwater lowland flooding',
    },
    'chalakudy': {
        'bbox': (76.28, 10.25, 76.38, 10.35),
        'category': 'flood',
        'label_type': 'sentinel1_sar',
        'pre_event': ('2018-07-01', '2018-07-25'),
        'post_event': ('2018-08-15', '2018-08-25'),
        'description': '2018 Chalakudy flood (Chalakkudy river basin)',
    },
    'thrissur_flood': {
        'bbox': (76.17, 10.47, 76.27, 10.57),
        'category': 'flood',
        'label_type': 'sentinel1_sar',
        'pre_event': ('2018-07-01', '2018-07-25'),
        'post_event': ('2018-08-15', '2018-08-25'),
        'description': '2018 Thrissur lowlands flood',
    },
    'pathanamthitta': {
        'bbox': (76.73, 9.22, 76.83, 9.32),
        'category': 'flood',
        'label_type': 'sentinel1_sar',
        'pre_event': ('2018-07-01', '2018-07-25'),
        'post_event': ('2018-08-15', '2018-08-25'),
        'description': '2018 Pathanamthitta flood',
    },
    'north_paravur': {
        'bbox': (76.17, 10.10, 76.27, 10.20),
        'category': 'flood',
        'label_type': 'sentinel1_sar',
        'pre_event': ('2018-07-01', '2018-07-25'),
        'post_event': ('2018-08-15', '2018-08-25'),
        'description': '2018 North Paravur flood',
    },
    'ranni': {
        'bbox': (76.73, 9.33, 76.83, 9.43),
        'category': 'flood',
        'label_type': 'sentinel1_sar',
        'pre_event': ('2018-07-01', '2018-07-25'),
        'post_event': ('2018-08-15', '2018-08-25'),
        'description': '2018 Ranni flood (Pampa river basin)',
    },
    'thiruvalla': {
        'bbox': (76.52, 9.33, 76.62, 9.43),
        'category': 'flood',
        'label_type': 'sentinel1_sar',
        'pre_event': ('2018-07-01', '2018-07-25'),
        'post_event': ('2018-08-15', '2018-08-25'),
        'description': '2018 Thiruvalla flood',
    },

    # =========================================================================
    # 3. SAFE BASELINE SITES (20)
    # =========================================================================
    'kottayam_city': {
        'bbox': (76.47, 9.54, 76.57, 9.64),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat urban area, no disaster history (Kottayam city)',
    },
    'ernakulam_city': {
        'bbox': (76.25, 9.93, 76.35, 10.03),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat urban area, no disaster history (Ernakulam city)',
    },
    'palakkad_city': {
        'bbox': (76.60, 10.73, 76.70, 10.83),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat urban area, no disaster history (Palakkad city)',
    },
    'thrissur_city': {
        'bbox': (76.16, 10.48, 76.26, 10.58),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat urban area, no disaster history (Thrissur city)',
    },
    'kozhikode_city': {
        'bbox': (75.73, 11.20, 75.83, 11.30),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat urban area, no disaster history (Kozhikode city)',
    },
    'kollam_city': {
        'bbox': (76.55, 8.84, 76.65, 8.94),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat urban area, no disaster history (Kollam city)',
    },
    'thiruvananthapuram': {
        'bbox': (76.89, 8.47, 76.99, 8.57),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat urban area, no disaster history (Thiruvananthapuram)',
    },
    'alappuzha_city': {
        'bbox': (76.29, 9.44, 76.39, 9.54),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat urban area, no disaster history (Alappuzha city)',
    },
    'kannur_city': {
        'bbox': (75.32, 11.82, 75.42, 11.92),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat urban area, no disaster history (Kannur city)',
    },
    'malappuram_city': {
        'bbox': (76.02, 10.99, 76.12, 11.09),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat urban area, no disaster history (Malappuram city)',
    },
    'kasaragod_town': {
        'bbox': (74.95, 12.45, 75.05, 12.55),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat urban area, no disaster history (Kasaragod town)',
    },
    'thalassery': {
        'bbox': (75.44, 11.70, 75.54, 11.80),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat coastal/urban area, no disaster history (Thalassery)',
    },
    'guruvayur': {
        'bbox': (75.99, 10.54, 76.09, 10.64),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat urban area, no disaster history (Guruvayur)',
    },
    'perinthalmanna': {
        'bbox': (76.18, 10.92, 76.28, 11.02),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Stable inland town, no disaster history (Perinthalmanna)',
    },
    'kayamkulam': {
        'bbox': (76.45, 9.12, 76.55, 9.22),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat coastal plain, no disaster history (Kayamkulam)',
    },
    'changanassery': {
        'bbox': (76.49, 9.39, 76.59, 9.49),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat urban area, no disaster history (Changanassery)',
    },
    'mavelikara': {
        'bbox': (76.50, 9.21, 76.60, 9.31),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat inland town, no disaster history (Mavelikara)',
    },
    'ottapalam': {
        'bbox': (76.33, 10.72, 76.43, 10.82),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Stable river valley town, no disaster history (Ottapalam)',
    },
    'kodungallur': {
        'bbox': (76.15, 10.18, 76.25, 10.28),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat coastal town, no disaster history (Kodungallur)',
    },
    'attingal': {
        'bbox': (76.77, 8.65, 76.87, 8.75),
        'category': 'safe',
        'label_type': 'none',
        'pre_event': None,
        'post_event': None,
        'description': 'Flat urban area, no disaster history (Attingal)',
    },
}


def get_sites_by_category(category: str) -> Dict[str, Dict[str, Any]]:
    """Retrieve all sites belonging to a specific category."""
    return {
        name: site for name, site in TRAINING_SITES.items()
        if site.get('category', '').lower() == category.lower()
    }


def get_all_sites() -> Dict[str, Dict[str, Any]]:
    """Retrieve the full dictionary of all 40 training sites."""
    return TRAINING_SITES


def get_site_bbox(site_name: str) -> Tuple[float, float, float, float]:
    """Retrieve the bounding box tuple for a specific site."""
    if site_name not in TRAINING_SITES:
        raise KeyError(f"Site '{site_name}' not found in TRAINING_SITES database.")
    return TRAINING_SITES[site_name]['bbox']
