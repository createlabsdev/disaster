// Constants & Configuration
const CONFIG = {
    DEFAULT_CENTER: [10.5, 76.2], // Kerala center
    DEFAULT_ZOOM: 7,
    API_BASE: '/api',
    REFRESH_RATES: {
        DAMS: 60 * 60 * 1000, // 60 mins
        ALERTS: 10 * 60 * 1000 // 10 mins
    },
    COLORS: {
        RED: '#ff4757',
        ORANGE: '#ffa502',
        BLUE: '#1e90ff',
        GREEN: '#2ed573'
    }
};

// Global State
let map;
let riskOverlayLayer = null;
let markersLayer = null;

// Initialize when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
    initMap();
    initSettings();
    initSearch();
    initLandingScreen();
    
    // Add map type dropdown listener
    const mapTypeSelect = document.getElementById('map-type-select');
    if (mapTypeSelect) {
        mapTypeSelect.addEventListener('change', (e) => {
            const selectedType = e.target.value;
            if (activeRiskLayer) map.removeLayer(activeRiskLayer);
            if (terrainRiskLayer) map.removeLayer(terrainRiskLayer);
            if (siltationRiskLayer) map.removeLayer(siltationRiskLayer);

            if (selectedType === 'active' && activeRiskLayer) {
                activeRiskLayer.addTo(map);
            } else if (selectedType === 'terrain' && terrainRiskLayer) {
                terrainRiskLayer.addTo(map);
            } else if (selectedType === 'siltation' && siltationRiskLayer) {
                siltationRiskLayer.addTo(map);
            }
            updateLegend(selectedType);
        });
    }
    
    // Initial data fetch
    fetchDams();
    fetchAlerts();
    
    // Setup intervals
    setInterval(fetchDams, CONFIG.REFRESH_RATES.DAMS);
    setInterval(fetchAlerts, CONFIG.REFRESH_RATES.ALERTS);
});

// Landing Page & Language Toggle Switch Logic
function initLandingScreen() {
    const landingScreen = document.getElementById('landing-screen');
    const enterBtn = document.getElementById('enter-dashboard-btn');
    const openGuideBtn = document.getElementById('open-guide-btn');
    const langOpts = document.querySelectorAll('.lang-opt');
    const langLabel = document.getElementById('lang-label-text');
    const mainTitle = document.getElementById('landing-main-title');
    const subtitle = document.getElementById('landing-subtitle');
    const ctaText = document.getElementById('cta-text');

    if (enterBtn && landingScreen) {
        enterBtn.addEventListener('click', () => {
            landingScreen.classList.add('hidden');
            if (map) {
                setTimeout(() => map.invalidateSize(), 350);
            }
        });
    }

    if (openGuideBtn && landingScreen) {
        openGuideBtn.addEventListener('click', () => {
            landingScreen.classList.remove('hidden');
        });
    }

    if (langOpts.length > 0) {
        langOpts.forEach(btn => {
            btn.addEventListener('click', () => {
                const lang = btn.getAttribute('data-lang');
                langOpts.forEach(b => b.classList.remove('active'));
                btn.classList.add('active');
                switchLanguage(lang);
            });
        });
    }

    function switchLanguage(lang) {
        const faqsQ = document.querySelectorAll('.faq-q');
        const faqsA = document.querySelectorAll('.faq-a');

        if (lang === 'ml') {
            if (mainTitle) mainTitle.textContent = "വർഷകാല സുരക്ഷാ മാർഗ്ഗനിർദ്ദേശങ്ങൾ";
            if (subtitle) subtitle.textContent = "കാലവർഷക്കെടുതികളിൽ പൗരന്മാർ അറിഞ്ഞിരിക്കേണ്ട അവശ്യ വിവരങ്ങൾ";
            if (ctaText) ctaText.textContent = "ലൈവ് ഡാഷ്‌ബോർഡിലേക്ക് പ്രവേശിക്കുക →";
            if (langLabel) langLabel.textContent = "ഭാഷ:";
        } else {
            if (mainTitle) mainTitle.textContent = "Monsoon Safety & Disaster Intelligence";
            if (subtitle) subtitle.textContent = "Essential answers for citizens during extreme rainfall & monsoon season";
            if (ctaText) ctaText.textContent = "Proceed to Live Dashboard →";
            if (langLabel) langLabel.textContent = "Language:";
        }

        faqsQ.forEach(el => {
            const text = el.getAttribute(`data-${lang}`);
            if (text) el.textContent = text;
        });

        faqsA.forEach(el => {
            const text = el.getAttribute(`data-${lang}`);
            if (text) el.textContent = text;
        });
    }
}

// Map Setup
function initMap() {
    map = L.map('map', {
        zoomControl: false // We'll add it in custom position
    }).setView(CONFIG.DEFAULT_CENTER, CONFIG.DEFAULT_ZOOM);

    // Add controls
    L.control.zoom({ position: 'bottomright' }).addTo(map);
    L.control.scale({ position: 'bottomleft', metric: true, imperial: false }).addTo(map);
    
    // Custom Locate Control
    const LocateControl = L.Control.extend({
        options: { position: 'bottomright' },
        onAdd: function() {
            const container = L.DomUtil.create('div', 'leaflet-bar leaflet-control leaflet-control-locate');
            const a = L.DomUtil.create('a', '', container);
            a.innerHTML = '📍';
            a.href = '#';
            a.title = 'Locate me';
            a.onclick = function(e) {
                e.preventDefault();
                map.locate({setView: true, maxZoom: 14});
            }
            return container;
        }
    });
    map.addControl(new LocateControl());

    // Tile Layers
    const darkLayer = L.tileLayer('https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png', {
        attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> &copy; <a href="https://carto.com/">CARTO</a>',
        subdomains: 'abcd',
        maxZoom: 20
    });

    const satelliteLayer = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}', {
        attribution: 'Tiles &copy; Esri &mdash; Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community',
        maxZoom: 19
    });
    
    const streetLayer = L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
        attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a>',
        maxZoom: 19
    });

    darkLayer.addTo(map);

    const baseMaps = {
        "Dark Theme": darkLayer,
        "Satellite": satelliteLayer,
        "Street Map": streetLayer
    };

    window.layerControl = L.control.layers(baseMaps, null, { position: 'topright' }).addTo(map);
    
    markersLayer = L.layerGroup().addTo(map);

    // Add Legend
    const legend = L.control({ position: 'bottomright' });
    legend.onAdd = function (map) {
        const div = L.DomUtil.create('div', 'info legend');
        div.innerHTML = `
            <h4>Risk Level</h4>
            <div class="legend-item"><span style="background: rgba(46, 204, 113, 0.9);"></span> Very Low</div>
            <div class="legend-item"><span style="background: rgba(168, 224, 108, 0.9);"></span> Low</div>
            <div class="legend-item"><span style="background: rgba(241, 196, 15, 0.9);"></span> Moderate</div>
            <div class="legend-item"><span style="background: rgba(230, 126, 34, 0.9);"></span> High</div>
            <div class="legend-item"><span style="background: rgba(231, 76, 60, 0.9);"></span> Very High</div>
        `;
        return div;
    };
    legend.addTo(map);
}

// Settings & LocalStorage
function getApiKey() {
    return localStorage.getItem('owm_api_key') || '';
}

function saveApiKey(key) {
    localStorage.setItem('owm_api_key', key.trim());
}

function initSettings() {
    const modal = document.getElementById('settings-modal');
    const btnOpen = document.getElementById('settings-btn');
    const btnClose = document.getElementById('close-modal-btn');
    const btnSave = document.getElementById('save-settings-btn');
    const inputKey = document.getElementById('owm-api-key');

    inputKey.value = getApiKey();

    btnOpen.addEventListener('click', () => modal.classList.remove('hidden'));
    btnClose.addEventListener('click', () => modal.classList.add('hidden'));
    
    modal.addEventListener('click', (e) => {
        if (e.target === modal) modal.classList.add('hidden');
    });

    btnSave.addEventListener('click', () => {
        const val = inputKey.value;
        if (val) {
            saveApiKey(val);
            modal.classList.add('hidden');
            showToast('Settings saved successfully', 'success');
        } else {
            showToast('API Key cannot be empty', 'error');
        }
    });
}

// Search & Predict Flow
function initSearch() {
    const form = document.getElementById('search-form');
    form.addEventListener('submit', async (e) => {
        e.preventDefault();
        const placeInput = document.getElementById('place-input').value.trim();
        const forecastSelect = document.getElementById('forecast-select');
        const forecastHours = forecastSelect ? parseInt(forecastSelect.value, 10) : 0;
        
        if (!placeInput) return;
        
        const apiKey = getApiKey();

        try {
            showLoading('Geocoding location...');
            
            // 1. Geocode
            const geocodeUrl = `${CONFIG.API_BASE}/geocode?place_name=${encodeURIComponent(placeInput)}`;
            const geoRes = await fetch(geocodeUrl);
            const geoData = await geoRes.json();
            
            if (!geoRes.ok) {
                throw new Error(geoData.detail || 'Location not found in Kerala. Please try a different name.');
            }
            
            const lat = parseFloat(geoData.lat);
            const lon = parseFloat(geoData.lon);
            
            // 2. Compute 10km bbox (approx 0.05 degrees)
            const west = lon - 0.05;
            const south = lat - 0.05;
            const east = lon + 0.05;
            const north = lat + 0.05;
            
            // Fly map
            map.flyTo([lat, lon], 13, { duration: 1.5 });
            
            // Wait for fly animation to mostly finish
            setTimeout(async () => {
                try {
                    showLoading('Running AI prediction pipeline... (new sites take 1-2 min)');
                    
                    // 3. API Request to backend
                    const response = await fetch(`${CONFIG.API_BASE}/predict`, {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({
                            place_name: placeInput,
                            west, south, east, north,
                            owm_api_key: apiKey,
                            forecast_hours: forecastHours
                        })
                    });
                    
                    if (!response.ok) {
                        let errMsg = 'Failed to analyze risk';
                        try {
                            const errData = await response.json();
                            errMsg = errData.detail || errMsg;
                        } catch { 
                            errMsg = await response.text();
                        }
                        throw new Error(errMsg);
                    }
                    
                    const data = await response.json();
                    
                    // 4. Update UI
                    const t = Date.now();
                    updateRiskMap(
                        `${data.active_map_url}?t=${t}`, 
                        `${data.terrain_map_url}?t=${t}`, 
                        data.siltation_map_url ? `${data.siltation_map_url}?t=${t}` : null,
                        bounds=[[south, west], [north, east]]
                    );
                    updateWeather(data.weather);
                    updateRiskGauge(data.active_risk_score, data.terrain_vulnerability, data.confidence_score);
                    
                    hideLoading();
                    showToast(`Analysis complete for ${placeInput} (${data.processing_time_s}s)`, 'success');
                } catch (err) {
                    hideLoading();
                    showToast(err.message, 'error');
                    console.error(err);
                }
            }, 1500);
            
        } catch (err) {
            hideLoading();
            showToast(err.message, 'error');
            console.error(err);
        }
    });
}

let activeRiskLayer = null;
let terrainRiskLayer = null;
let siltationRiskLayer = null;

function updateRiskMap(activeUrl, terrainUrl, siltationUrl, bounds) {
    if (activeRiskLayer) map.removeLayer(activeRiskLayer);
    if (terrainRiskLayer) map.removeLayer(terrainRiskLayer);
    if (siltationRiskLayer) map.removeLayer(siltationRiskLayer);
    
    if (terrainUrl) {
        terrainRiskLayer = L.imageOverlay(terrainUrl, bounds, {
            opacity: 0.80,
            interactive: false
        });
    }

    if (activeUrl) {
        activeRiskLayer = L.imageOverlay(activeUrl, bounds, {
            opacity: 0.80,
            interactive: false
        });
    }

    if (siltationUrl) {
        siltationRiskLayer = L.imageOverlay(siltationUrl, bounds, {
            opacity: 0.85,
            interactive: false
        });
    }
    
    // Show the currently selected map type
    const selectEl = document.getElementById('map-type-select');
    const selectedType = selectEl ? selectEl.value : 'active';
    
    if (selectedType === 'active' && activeRiskLayer) {
        activeRiskLayer.addTo(map);
    } else if (selectedType === 'terrain' && terrainRiskLayer) {
        terrainRiskLayer.addTo(map);
    } else if (selectedType === 'siltation' && siltationRiskLayer) {
        siltationRiskLayer.addTo(map);
    }
    updateLegend(selectedType);
}

function updateLegend(mapType) {
    const legendDiv = document.querySelector('.info.legend');
    if (!legendDiv) return;

    if (mapType === 'siltation') {
        legendDiv.innerHTML = `
            <h4 style="margin-bottom:6px; color:#f1c40f;">Silt & Mud Deposit</h4>
            <div class="legend-item"><span style="background: rgba(241, 196, 15, 0.9);"></span> Low Siltation</div>
            <div class="legend-item"><span style="background: rgba(230, 126, 34, 0.9);"></span> Moderate Deposit (25-40%)</div>
            <div class="legend-item"><span style="background: rgba(192, 57, 43, 0.9);"></span> Severe Siltation (>45%)</div>
            <div style="margin-top:6px; font-size:0.72rem; color:#bbb; max-width:160px; line-height:1.2;">
                ⚠️ <i>High siltation reduces river carrying capacity & causes early overflow.</i>
            </div>
        `;
    } else {
        legendDiv.innerHTML = `
            <h4>Risk Level</h4>
            <div class="legend-item"><span style="background: rgba(46, 204, 113, 0.9);"></span> Very Low</div>
            <div class="legend-item"><span style="background: rgba(168, 224, 108, 0.9);"></span> Low</div>
            <div class="legend-item"><span style="background: rgba(241, 196, 15, 0.9);"></span> Moderate</div>
            <div class="legend-item"><span style="background: rgba(230, 126, 34, 0.9);"></span> High</div>
            <div class="legend-item"><span style="background: rgba(231, 76, 60, 0.9);"></span> Very High</div>
        `;
    }
}

// Weather Update
function updateWeather(data) {
    if (!data) return;
    
    document.getElementById('val-temp').textContent = data.temp_c != null ? `${Math.round(data.temp_c)}°C` : '--°C';
    document.getElementById('val-humidity').textContent = data.humidity != null ? `${Math.round(data.humidity)}%` : '--%';
    document.getElementById('val-rain').textContent = `${data.rainfall_mm_h || 0} mm/h`;
    document.getElementById('val-wind').textContent = data.wind_speed_kmh != null ? `${data.wind_speed_kmh} km/h` : '-- km/h';
    
    const descEl = document.getElementById('weather-desc');
    let emoji = '🌤️';
    const main = (data.description || '').toLowerCase();
    
    if (main.includes('storm') || main.includes('thunder')) emoji = '⛈️';
    else if (main.includes('rain') || main.includes('drizzle')) emoji = '🌧️';
    else if (main.includes('snow')) emoji = '❄️';
    else if (main.includes('fog')) emoji = '🌫️';
    else if (main.includes('cloud') || main.includes('overcast')) emoji = '☁️';
    else if (main.includes('clear')) emoji = '☀️';
    
    const dailyTotalStr = data.daily_total_mm !== undefined ? ` | Daily Total: ${data.daily_total_mm} mm` : '';
    descEl.textContent = `${emoji} ${data.description ? data.description.charAt(0).toUpperCase() + data.description.slice(1) : 'Unknown'}${dailyTotalStr}`;
}

// Risk Gauge Animate
function updateRiskGauge(activeScore, terrainScore, confidence) {
    activeScore = Math.max(0, Math.min(100, activeScore || 0));
    terrainScore = Math.max(0, Math.min(100, terrainScore || 0));
    const confVal = (confidence != null) ? Math.round(confidence) : '--';
    
    // Active Risk UI Elements
    const activeFill = document.getElementById('risk-gauge-fill');
    const activeValText = document.getElementById('risk-score-value');
    const activeLevelText = document.getElementById('risk-level-text');
    
    // Terrain Vulnerability UI Elements
    const terrainFill = document.getElementById('terrain-gauge-fill');
    const terrainValText = document.getElementById('terrain-score-value');
    
    const confText = document.getElementById('confidence-score-value');
    if (confText) confText.textContent = confVal;
    
    // Circle circumference is 314 (2 * pi * r = 2 * 3.14 * 50)
    const circumference = 314;
    
    // 1. Update Active Risk Gauge
    const activeOffset = circumference - (activeScore / 100) * circumference;
    activeFill.style.strokeDashoffset = activeOffset;
    
    let activeCurrent = parseInt(activeValText.textContent) || 0;
    const activeStep = activeScore > activeCurrent ? 1 : -1;
    const activeTimer = setInterval(() => {
        if (activeCurrent === Math.round(activeScore)) {
            clearInterval(activeTimer);
        } else {
            activeCurrent += activeStep;
            activeValText.textContent = activeCurrent;
        }
    }, 20);
    
    let activeColor = CONFIG.COLORS.GREEN;
    let activeText = 'LOW RISK';
    
    if (activeScore >= 70) {
        activeColor = CONFIG.COLORS.RED;
        activeText = 'CRITICAL RISK';
    } else if (activeScore >= 30) {
        activeColor = CONFIG.COLORS.ORANGE;
        activeText = 'MODERATE RISK';
    }
    
    activeFill.style.stroke = activeColor;
    activeLevelText.textContent = activeText;
    activeLevelText.style.color = activeColor;

    // 2. Update Terrain Vulnerability Gauge
    const terrainOffset = circumference - (terrainScore / 100) * circumference;
    terrainFill.style.strokeDashoffset = terrainOffset;
    
    let terrainCurrent = parseInt(terrainValText.textContent) || 0;
    const terrainStep = terrainScore > terrainCurrent ? 1 : -1;
    const terrainTimer = setInterval(() => {
        if (terrainCurrent === Math.round(terrainScore)) {
            clearInterval(terrainTimer);
        } else {
            terrainCurrent += terrainStep;
            terrainValText.textContent = terrainCurrent;
        }
    }, 20);
}

// Dams API Fetch & Render
async function fetchDams() {
    try {
        const cached = localStorage.getItem('dams_data');
        const cachedTime = localStorage.getItem('dams_time');
        const now = Date.now();
        
        if (cached && cachedTime && (now - parseInt(cachedTime) < CONFIG.REFRESH_RATES.DAMS)) {
            renderDams(JSON.parse(cached));
            // Fetch in background quietly
            fetchDamsQuiet();
            return;
        }
        
        await fetchDamsQuiet();
    } catch (err) {
        console.error("Dam fetch error:", err);
        document.getElementById('dams-list').innerHTML = '<div class="loading-text" style="color:var(--accent-red)">Failed to load dam data</div>';
    }
}

async function fetchDamsQuiet() {
    try {
        const res = await fetch(`${CONFIG.API_BASE}/dams`);
        if (!res.ok) throw new Error('API error');
        const data = await res.json();
        
        localStorage.setItem('dams_data', JSON.stringify(data));
        localStorage.setItem('dams_time', Date.now().toString());
        
        renderDams(data);
    } catch (err) {
        console.warn("Quiet dam fetch failed", err);
    }
}

function renderDams(dams) {
    const list = document.getElementById('dams-list');
    list.innerHTML = '';
    markersLayer.clearLayers();
    
    if (!dams || dams.length === 0) {
        list.innerHTML = '<div class="loading-text">No dam data available</div>';
        return;
    }
    
    // Sort by alert severity: RED, ORANGE, BLUE, NORMAL (GREEN)
    const severityOrder = { 'RED': 0, 'ORANGE': 1, 'BLUE': 2, 'NORMAL': 3 };
    dams.sort((a, b) => {
        const sA = severityOrder[(a.alert || 'NORMAL').toUpperCase()] ?? 4;
        const sB = severityOrder[(b.alert || 'NORMAL').toUpperCase()] ?? 4;
        return sA - sB;
    });
    
    dams.forEach(dam => {
        const alertLvl = (dam.alert || 'NORMAL').toUpperCase();
        let color = CONFIG.COLORS.GREEN;
        if (alertLvl === 'RED') color = CONFIG.COLORS.RED;
        else if (alertLvl === 'ORANGE') color = CONFIG.COLORS.ORANGE;
        else if (alertLvl === 'BLUE') color = CONFIG.COLORS.BLUE;
        
        // HTML Card
        const card = document.createElement('div');
        card.className = 'dam-card';
        card.style.borderLeftColor = color;
        
        card.innerHTML = `
            <div class="dam-info">
                <h3>${dam.name || 'Unknown Dam'}</h3>
                <p>${dam.district || 'Kerala'}</p>
            </div>
            <div class="dam-level">
                ${formatNumber(dam.current_level)} m
                <span>Full: ${formatNumber(dam.full_level)} m</span>
                <div class="badge" style="background:${color}33; color:${color}">${alertLvl}</div>
            </div>
        `;
        
        list.appendChild(card);
        
        // Map Marker
        if (dam.lat && dam.lon) {
            L.circleMarker([dam.lat, dam.lon], {
                radius: 8,
                fillColor: color,
                color: '#fff',
                weight: 1,
                opacity: 1,
                fillOpacity: 0.8
            }).bindPopup(`<b>${dam.name}</b><br>Level: ${dam.current_level}m / ${dam.full_level}m<br>Status: ${alertLvl}`)
              .addTo(markersLayer);
        }
    });
}

// Alerts API Fetch & Render
async function fetchAlerts() {
    try {
        const cached = localStorage.getItem('alerts_data');
        const cachedTime = localStorage.getItem('alerts_time');
        const now = Date.now();
        
        if (cached && cachedTime && (now - parseInt(cachedTime) < CONFIG.REFRESH_RATES.ALERTS)) {
            renderAlerts(JSON.parse(cached));
            fetchAlertsQuiet();
            return;
        }
        await fetchAlertsQuiet();
    } catch (err) {
        console.error("Alert fetch error:", err);
        document.getElementById('alerts-list').innerHTML = '<div class="loading-text" style="color:var(--accent-red)">Failed to load alerts</div>';
    }
}

async function fetchAlertsQuiet() {
    try {
        const res = await fetch(`${CONFIG.API_BASE}/alerts`);
        if (!res.ok) throw new Error('API error');
        const data = await res.json();
        
        localStorage.setItem('alerts_data', JSON.stringify(data));
        localStorage.setItem('alerts_time', Date.now().toString());
        
        renderAlerts(data);
    } catch (err) {
        console.warn("Quiet alert fetch failed", err);
    }
}

function renderAlerts(alerts) {
    const list = document.getElementById('alerts-list');
    list.innerHTML = '';
    
    if (!alerts || alerts.length === 0) {
        list.innerHTML = '<div class="loading-text">No active alerts from IMD.</div>';
        return;
    }
    
    alerts.forEach((alert, index) => {
        const severity = (alert.severity || '').toLowerCase();
        let sClass = '';
        let icon = 'ℹ️';
        let color = CONFIG.COLORS.BLUE;
        
        if (severity === 'red') { sClass = 'severity-extreme'; icon = '🚨'; color = CONFIG.COLORS.RED; }
        else if (severity === 'orange') { icon = '⚠️'; color = CONFIG.COLORS.ORANGE; }
        else if (severity === 'yellow') { icon = '⚡'; color = CONFIG.COLORS.ORANGE; }
        
        const item = document.createElement('div');
        item.className = `alert-item ${sClass}`;
        
        item.innerHTML = `
            <div class="alert-header">
                <div class="alert-title">
                    <span style="color:${color}">${icon}</span>
                    ${alert.title || 'Alert'}
                </div>
                <svg class="alert-chevron" xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polyline points="6 9 12 15 18 9"></polyline></svg>
            </div>
            <div class="alert-body">
                <p><strong>Published:</strong> ${alert.published || 'N/A'}</p>
                <p style="margin-top:0.5rem">${alert.description || 'No description provided.'}</p>
            </div>
        `;
        
        // Accordion functionality
        item.querySelector('.alert-header').addEventListener('click', () => {
            const isActive = item.classList.contains('active');
            // Close others
            document.querySelectorAll('.alert-item').forEach(el => el.classList.remove('active'));
            if (!isActive) item.classList.add('active');
        });
        
        // Open first one by default
        if (index === 0) item.classList.add('active');
        
        list.appendChild(item);
    });
}

// Utilities
function showLoading(text) {
    const overlay = document.getElementById('loading-overlay');
    const textEl = document.getElementById('loading-text');
    textEl.textContent = text || 'Loading...';
    overlay.classList.remove('hidden');
}

function hideLoading() {
    const overlay = document.getElementById('loading-overlay');
    overlay.classList.add('hidden');
}

function showToast(message, type = 'info') {
    const container = document.getElementById('toast-container');
    const toast = document.createElement('div');
    toast.className = `toast ${type}`;
    toast.textContent = message;
    
    container.appendChild(toast);
    
    setTimeout(() => {
        toast.style.animation = 'slideInRight 0.3s ease reverse forwards';
        setTimeout(() => toast.remove(), 300);
    }, 4000);
}

function formatNumber(n) {
    if (n === null || n === undefined) return '--';
    return Number(n).toLocaleString('en-IN', { maximumFractionDigits: 2 });
}
