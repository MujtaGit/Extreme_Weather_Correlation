import cdsapi
import datetime
import pandas as pd

# Initialize the API client
c = cdsapi.Client()

def read_events_from_csv(file_path):
    return pd.read_csv(file_path)

# Function to retrieve data for a specific event
def retrieve_event_data(year, month, month_end, start_day, end_day, file_name):
    # Generate a date range spanning multiple months if necessary
    start_date = datetime.datetime(year, month, start_day)
    end_date = datetime.datetime(year, month_end if month_end else month, end_day)
    
    current_date = start_date
    while current_date <= end_date:
        days = []
        month = current_date.month
        while current_date <= end_date and current_date.month == month:
            days.append(current_date.strftime("%d"))
            current_date += datetime.timedelta(days=1)
        
        # Generate a unique file name for each month to avoid overwriting
        monthly_file_name = f"{file_name}_{str(month).zfill(2)}.grib"

        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'format': 'grib',
                'variable': [
                    '10m_u_component_of_wind', '10m_v_component_of_wind', 'total_precipitation',
                    'surface_pressure', 'sea_surface_temperature'
                ],
                'year': str(year),
                'month': str(month).zfill(2),
                'day': days,
                'time': [
                    '00:00', '01:00', '02:00', '03:00', '04:00', '05:00',
                    '06:00', '07:00', '08:00', '09:00', '10:00', '11:00',
                    '12:00', '13:00', '14:00', '15:00', '16:00', '17:00',
                    '18:00', '19:00', '20:00', '21:00', '22:00', '23:00',
                ],
                'area': [
                    22, -95, 32, -70,  # Adjust this based on the specific area of interest
                ]
            },
            monthly_file_name
        )



# Loop through each event to retrieve the data
# Loop and retrieve data

def main():  
    file_path = "Events_temp.csv"  
    events = read_events_from_csv(file_path)  
  
    # Loop through each event to retrieve the data  
    for index, event in events.iterrows():
        file_name = f"data_temp/{event['name'].replace(' ', '_')}"
        retrieve_event_data(event['year'], event['month'], event.get('month_end'), event['start_day'], event['end_day'], file_name)
        
if __name__ == "__main__":  
    main()