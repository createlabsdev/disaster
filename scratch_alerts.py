import httpx
import asyncio
import xml.etree.ElementTree as ET

async def test():
    urls = [
        "https://sachet.ndma.gov.in/cap_public_website/rsspage.php",
        "https://sachet.ndma.gov.in/cap_public_website/rss_state.php?state_code=KL"
    ]
    for url in urls:
        print(f"Fetching {url}...")
        try:
            async with httpx.AsyncClient() as client:
                resp = await client.get(url, timeout=10, follow_redirects=True)
                print(f"Status: {resp.status_code}")
                if resp.status_code == 200:
                    root = ET.fromstring(resp.text)
                    channel = root.find("channel")
                    items = channel.findall("item")
                    print(f"Found {len(items)} items")
                    for i, item in enumerate(items[:5]):
                        t = item.findtext("title", "")
                        d = item.findtext("description", "")
                        print(f" - [{i}] {t}")
                        print(f"   Desc: {d[:100]}")
        except Exception as e:
            print(f"Error: {e}")

asyncio.run(test())
