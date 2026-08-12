import requests
import sys

def geocode(place):
    url = "https://nominatim.openstreetmap.org/search"
    params = {
        "q": place,
        "format": "json",
        "limit": 1,
        "countrycodes": "in"
    }
    headers = {"User-Agent": "KeralaDisasterDashboard/1.0"}
    resp = requests.get(url, params=params, headers=headers)
    data = resp.json()
    if not data:
        print(f"'{place}' -> Not found")
        return
    print(f"'{place}' -> {data[0]['display_name']}")

geocode("Delhi")
geocode("Munnar")
geocode("Kottayam")
geocode("New Delhi")
geocode("Ernakulam")
